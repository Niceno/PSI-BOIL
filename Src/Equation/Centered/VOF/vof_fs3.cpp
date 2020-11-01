#include "vof.h"

/******************************************************************************/
void VOF::cal_fs3(const Scalar & scp) {
/***************************************************************************//**
 \brief Calculate free-surface position between cell centers
    if there is no interface in the cell, unreal=yotta (=1e+24) is stored.
    plane: vm1*x + vm2*y + vm3*z = alpha
    output: fs
*******************************************************************************/

  /* tolerance is necessary because the linear-planes are not closed */
  /* with this version, it shouldnt be necessary anymore */
  real tolf = 0.0e-2;

  /* fs is a result of calculations involving nx,ny,nz,nalpha,color
   * and cell dimensions. Thus, by taking care of properly exchanging
   * these variables, we can safely venture to calculate fs in buffer
   * cells without dealing with the exchange() method, which would 
   * not function properly for fs. So don't worry and let's loop over
   * (almost) all cells, including buffers except for the furthest one */
  int ibeg(si()), iend(ei()+1);
  int jbeg(sj()), jend(ej()+1);
  int kbeg(sk()), kend(ek()+1);

  if(bflag_struct.ifull) {
    if(bflag_struct.iminp) {
      ibeg -= boil::BW-1;
    }
    if(bflag_struct.imaxp) {
      iend += boil::BW-1;
    }
  }
  if(bflag_struct.jfull) {
    if(bflag_struct.jminp) {
      jbeg -= boil::BW-1;
    }
    if(bflag_struct.jmaxp) {
      jend += boil::BW-1;
    }
  }
  if(bflag_struct.kfull) {
    if(bflag_struct.kminp) {
      kbeg -= boil::BW-1;
    }
    if(bflag_struct.kmaxp) {
      kend += boil::BW-1;
    }
  }
        
  /* initialize */
  for_m(m)
    for_avmijk(fs,m,i,j,k)
      fs[m][i][j][k] = boil::unreal;

  Comp m;

  /******************************************
  *             x-direction                 *
  ******************************************/

  m = Comp::i();
  for(int i=ibeg; i<=iend; i++)
  for(int j=sj(); j<=ej(); j++)
  for(int k=sk(); k<=ek(); k++) {

    /* degenerate cases */
    real clrw = scp[i-1][j][k];
    real clre = scp[i  ][j][k];

    if((clrw-phisurf)*(clre-phisurf)>0.0) 
      continue;
  
    if(  (clrw<boil::pico&&clre-1.0>-boil::pico)
       ||(clrw-1.0>-boil::pico&&clre<boil::pico)) {
      fs[m][i][j][k] = scp.xn(i);
      continue;
    }

    /* calculate normalized candidate positions */
    real fsxw = fs_val(m,i-1,j,k);
    real fsxe = fs_val(m,i  ,j,k);

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flagw = (0.5     <= fsxw && fsxw <= 1.0+tolf);
    bool flage = (0.0-tolf <= fsxe && fsxe <= 0.5    );

    if     ( flagw && !flage) { /* west is real and east is not */
      fs[m][i][j][k] = scp.xn(i-1) + scp.dxc(i-1) * fsxw;
      continue;
    }
    else if(!flagw &&  flage) { /* east is real and west is not */
      fs[m][i][j][k] = scp.xn(i) + scp.dxc(i) * fsxe;
      continue;
    }
    else if( flagw &&  flage) { /* both are real */

      /* calculate distance between the two candidate positions */
      real fsx_diff = fsxe + 1.0 - fsxw;
      if(fabs(fsx_diff)<tolf) { /* they are close to each other */
        fs[m][i][j][k] = scp.xn(i);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = (1.0-fsxw>fsxe) ? 
                         scp.xn(i-1) + scp.dxc(i-1) * fsxw
                       : scp.xn(i  ) + scp.dxc(i  ) * fsxe;
      }
      continue;
    }

    fs[m][i][j][k] = scp.xn(i);
  } 
  
  /******************************************
  *             y-direction                 *
  ******************************************/

  m = Comp::j();
  for(int i=si(); i<=ei(); i++)
  for(int j=jbeg; j<=jend; j++)
  for(int k=sk(); k<=ek(); k++) {

    /* degenerate cases */
    real clrs = scp[i][j-1][k];
    real clrn = scp[i][j  ][k];

    if((clrs-phisurf)*(clrn-phisurf)>0.0) 
      continue;
    
    if(  (clrs<boil::pico&&clrn-1.0>-boil::pico)
       ||(clrs-1.0>-boil::pico&&clrn<boil::pico)) {
      fs[m][i][j][k] = scp.yn(j);
      continue;
    }

    /* calculate normalized candidate positions */
    real fsys = fs_val(m,i,j-1,k);
    real fsyn = fs_val(m,i,j  ,k);

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flags = (0.5     <= fsys && fsys <= 1.0+tolf);
    bool flagn = (0.0-tolf <= fsyn && fsyn <= 0.5);
   
    if     ( flags && !flagn) { /* south is real and north is not */
      fs[m][i][j][k] = scp.yn(j-1) + scp.dyc(j-1) * fsys;
      continue;
    }
    else if(!flags &&  flagn) { /* north is real and south is not */
      fs[m][i][j][k] = scp.yn(j) + scp.dyc(j) * fsyn;
      continue;
    }
    else if( flags &&  flagn) { /* both are real */

      /* calculate distance between the two candidate positions */
      real fsy_diff = fsyn + 1.0 - fsys;
      if(fabs(fsy_diff)<tolf) { /* they are close to each other */
        fs[m][i][j][k] = scp.yn(j);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = (1.0-fsys>fsyn) ? 
                         scp.yn(j-1) + scp.dyc(j-1) * fsys
                       : scp.yn(j  ) + scp.dyc(j  ) * fsyn;
      }
      continue;
    }
    fs[m][i][j][k] = scp.yn(j);
  } 
  
  /******************************************
  *             z-direction                 *
  ******************************************/

  m = Comp::k();
  for(int i=si(); i<=ei(); i++)
  for(int j=sj(); j<=ej(); j++)
  for(int k=kbeg; k<=kend; k++) {

    /* degenerate cases */
    real clrb = scp[i][j][k-1];
    real clrt = scp[i][j][k  ];

    if((clrb-phisurf)*(clrt-phisurf)>0.0) 
      continue;
    
    if(  (clrb<boil::pico&&clrt-1.0>-boil::pico)
       ||(clrb-1.0>-boil::pico&&clrt<boil::pico)) {
      fs[m][i][j][k] = scp.zn(k);
      continue;
    }

    /* calculate normalized candidate positions */
    real fszb = fs_val(m,i,j,k-1);
    real fszt = fs_val(m,i,j,k  );

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flagb = (0.5     <= fszb && fszb <= 1.0+tolf);
    bool flagt = (0.0-tolf <= fszt && fszt <= 0.5    );
   
    //if(i==46&&j==50&&k==76) boil::oout<<"VOF::cal_FS3 "<<fszb<<" "<<fszt<<boil::endl;
  
    if     ( flagb && !flagt) { /* bottom is real and top is not */
      fs[m][i][j][k] = scp.zn(k-1) + scp.dzc(k-1) * fszb;
      continue;
    }
    else if(!flagb &&  flagt) { /* top is real and bottom is not */
      fs[m][i][j][k] = scp.zn(k) + scp.dzc(k) * fszt;
      continue;
    }
    else if( flagb &&  flagt) { /* both are real */

      /* calculate distance between the two candidate positions */
      real fsz_diff = fszt + 1.0 - fszb;
      if(fabs(fsz_diff)<tolf) { /* they are close to each other */
        fs[m][i][j][k] = scp.zn(k);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = (1.0-fszb>fszt) ? 
                         scp.zn(k-1) + scp.dzc(k-1) * fszb
                       : scp.zn(k  ) + scp.dzc(k  ) * fszt;
      }
      continue;
    }
    fs[m][i][j][k] = scp.zn(k);
  } 

  /* correct at boundaries */
  if       (subgrid_method == SubgridMethod::PLIC()) {
    heavi->fs_bnd_geometric(scp,fs,tol_wall);
  } else if(subgrid_method == SubgridMethod::None()) {
    heavi->fs_bnd_nosubgrid(scp,fs,tol_wall);
  } else if(subgrid_method == SubgridMethod::SLICliquid()) {
    heavi->fs_bnd_1D(scp,fs,tol_wall,Sign::pos());
  } else if(subgrid_method == SubgridMethod::SLICgas()) {
    heavi->fs_bnd_1D(scp,fs,tol_wall,Sign::neg());
  } else {
    boil::oout<<"VOF::fs3: Subgrid method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }
  /* symmetry */
  fs_bnd_symmetry(scp,fs,tol_wall);
  //fs.exchange_all();

  //boil::plot->plot(fs,scp, "fs-clr", 0);
 
  return;
}

/***********************
 * ancillary function
 ***********************/
real VOF::fs_val(const Comp m, const int i, const int j, const int k) {

  real alpha = nalpha[i][j][k];
  /* degenerate case */
  if(!boil::realistic(alpha))
    return boil::unreal;

  real vn1 = -nx[i][j][k];
  real vn2 = -ny[i][j][k];
  real vn3 = -nz[i][j][k];

  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real xpos = 0.5;
  real ypos = 0.5;
  real zpos = 0.5;

#if 1
  if(m==Comp::i()) {
    real xuni = (alpha-vm2*ypos-vm3*zpos)/(vm1+boil::pico);
    if(vn1<0)
      xuni = 1.0-xuni;
    return xuni;
  } else if(m==Comp::j()) {
    real yuni = (alpha-vm1*xpos-vm3*zpos)/(vm2+boil::pico);
    if(vn2<0)
      yuni = 1.0-yuni;
    return yuni;
  } else {
    real zuni = (alpha-vm1*xpos-vm2*ypos)/(vm3+boil::pico);
    if(vn3<0)
      zuni = 1.0-zuni;
    return zuni;
  }
#else
  if(m==Comp::i()) {
    real xuni = (alpha-vn2*ypos-vn3*zpos)/(vm1+boil::pico);
    if(vn1<0) 
      xuni = xuni-1.0;
    return xuni;
  } else if(m==Comp::j()) {
    real yuni = (alpha-vn1*xpos-vn3*zpos)/(vm2+boil::pico);
    if(vn2<0)
      yuni = yuni-1.0;
    return yuni;
  } else {
    real zuni = (alpha-vn1*xpos-vn2*ypos)/(vm3+boil::pico);
    if(vn3<0)
      zuni = zuni-1.0;
    return zuni;
  }
#endif

  return boil::unreal;
}

