#include "heaviside.h"

/******************************************************************************/
void Heaviside::cal_fs_geom(const Scalar & scp, 
                            const Scalar & nx, const Scalar & ny,
                            const Scalar & nz, const Scalar & nalpha,
                            Vector & fs, 
                            const real tol_wall, const bool use_subgrid) {
/***************************************************************************//**
 \brief Calculate free-surface position between cell centers
    if there is no interface in the cell, unreal=yotta (=1e+24) is stored.
    plane: vm1*x + vm2*y + vm3*z = alpha
    output: fs
*******************************************************************************/

  /* tolerance, not really necessary */
  real tolf = 0.0e-2;
  /* if you decide to use it, don't forget to scale by grid spacing, as this is
     currently implemented as a constant tolerance! */

  /* initialize */
  for_m(m)
    for_avmijk(fs,m,i,j,k)
      fs[m][i][j][k] = boil::unreal;

  Comp m;
  
  /******************************************
  *             x-direction                 *
  ******************************************/

  m = Comp::i();
  for(int i=scp.si(); i<=scp.ei()+1; i++)
  for(int j=scp.sj(); j<=scp.ej()  ; j++)
  for(int k=scp.sk(); k<=scp.ek()  ; k++) {
  
    /* degenerate cases */
    real clrw = scp[i-1][j][k];
    real clre = scp[i  ][j][k];

    if((clrw-clrsurf)*(clre-clrsurf)>0.0) 
      continue;
  
    if(  (clrw<boil::pico&&clre-1.0>-boil::pico)
       ||(clrw-1.0>-boil::pico&&clre<boil::pico)) {
      fs[m][i][j][k] = scp.xn(i);
      boil::oout<<i<<" "<<j<<" "<<k<<" "<<clrw<<" "<<clre<<boil::endl;
      exit(0);
      continue;
    }

    /* calculate candidate positions */
    real fsxw = fs_val(m,i-1,j,k,nx,ny,nz,nalpha);
    real fsxe = fs_val(m,i  ,j,k,nx,ny,nz,nalpha);

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flagw = (scp.xc(i-1)    <= fsxw && fsxw <= scp.xn(i)+tolf);
    bool flage = (scp.xn(i)-tolf <= fsxe && fsxe <= scp.xc(i)     );

    if     ( flagw && !flage) { /* west is real and east is not */
      fs[m][i][j][k] = fsxw;
      continue;
    }
    else if(!flagw &&  flage) { /* east is real and west is not */
      fs[m][i][j][k] = fsxe;
      continue;
    }
    else if( flagw &&  flage) { /* both are real */

      /* calculate distance between the two candidate positions */
      real fsx_diff = fsxe - fsxw;
      if(fabs(fsx_diff)<tolf) { /* they are close to each other */
        fs[m][i][j][k] = scp.xn(i);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = fabs(fsxw-scp.xn(i))>fabs(fsxe-scp.xn(i)) 
                         ? fsxw : fsxe;
      }
      continue;
    }

      boil::oout<<i<<" "<<j<<" "<<k<<" "<<flagw<<" "<<flage<<" | "<<fsxw<<" "<<fsxe<<" | "<<scp.xn(i)<<" "<<scp.xc(i-1)<<" | "<<clrw<<" "<<clre<<boil::endl;
      //exit(0);
    fs[m][i][j][k] = scp.xn(i);
  } 
  
  /******************************************
  *             y-direction                 *
  ******************************************/

  m = Comp::j();
  for(int i=scp.si(); i<=scp.ei()  ; i++)
  for(int j=scp.sj(); j<=scp.ej()+1; j++)
  for(int k=scp.sk(); k<=scp.ek()  ; k++) {
               
    /* degenerate cases */
    real clrs = scp[i][j-1][k];
    real clrn = scp[i][j  ][k];

    if((clrs-clrsurf)*(clrn-clrsurf)>0.0) 
      continue;
    
    if(  (clrs<boil::pico&&clrn-1.0>-boil::pico)
       ||(clrs-1.0>-boil::pico&&clrn<boil::pico)) {
      fs[m][i][j][k] = scp.yn(j);
      continue;
    }

    /* calculate  candidate positions */
    real fsys = fs_val(m,i,j-1,k,nx,ny,nz,nalpha);
    real fsyn = fs_val(m,i,j  ,k,nx,ny,nz,nalpha);

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flags = (scp.yc(j-1)    <= fsys && fsys <= scp.yn(j)+tolf);
    bool flagn = (scp.yn(j)-tolf <= fsyn && fsyn <= scp.yc(j)     );
   
    if     ( flags && !flagn) { /* south is real and north is not */
      fs[m][i][j][k] = fsys;
      continue;
    }
    else if(!flags &&  flagn) { /* north is real and south is not */
      fs[m][i][j][k] = fsyn;
      continue;
    }
    else if( flags &&  flagn) { /* both are real */

      /* calculate distance between the two candidate positions */
      real fsy_diff = fsyn - fsys;
      if(fabs(fsy_diff)<tolf) { /* they are close to each other */
        fs[m][i][j][k] = scp.yn(j);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = fabs(fsys-scp.yn(j))>fabs(fsyn-scp.yn(j)) 
                         ? fsys : fsyn;
      }
      continue;
    }
    fs[m][i][j][k] = scp.yn(j);
  } 
  
  /******************************************
  *             z-direction                 *
  ******************************************/

  m = Comp::k();
  for(int i=scp.si(); i<=scp.ei()  ; i++)
  for(int j=scp.sj(); j<=scp.ej()  ; j++)
  for(int k=scp.sk(); k<=scp.ek()+1; k++) {

    /* degenerate cases */
    real clrb = scp[i][j][k-1];
    real clrt = scp[i][j][k  ];

    if((clrb-clrsurf)*(clrt-clrsurf)>0.0) 
      continue;
    
    if(  (clrb<boil::pico&&clrt-1.0>-boil::pico)
       ||(clrb-1.0>-boil::pico&&clrt<boil::pico)) {
      fs[m][i][j][k] = scp.zn(k);
      continue;
    }

    /* calculate candidate positions */
    real fszb = fs_val(m,i,j,k-1,nx,ny,nz,nalpha);
    real fszt = fs_val(m,i,j,k  ,nx,ny,nz,nalpha);

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flagb = (scp.zc(k-1)    <= fszb && fszb <= scp.zn(k)+tolf);
    bool flagt = (scp.zn(k)-tolf <= fszt && fszt <= scp.zc(k)     );
   
    if     ( flagb && !flagt) { /* bottom is real and top is not */
      fs[m][i][j][k] = fszb;
      continue;
    }
    else if(!flagb &&  flagt) { /* top is real and bottom is not */
      fs[m][i][j][k] = fszt;
      continue;
    }
    else if( flagb &&  flagt) { /* both are real */

      /* calculate distance between the two candidate positions */
      real fsz_diff = fszt - fszb;
      if(fabs(fsz_diff)<tolf) { /* they are close to each other */
        fs[m][i][j][k] = scp.zn(k);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = fabs(fszb-scp.zn(k))>fabs(fszt-scp.zn(k)) 
                         ? fszb : fszt;
      }
      continue;
    }
    fs[m][i][j][k] = scp.zn(k);
  } 

  /* correct at boundaries */
  if(use_subgrid)
    fs_bnd_subgrid(scp,fs,tol_wall);
  else
    fs_bnd_nosubgrid(scp,fs,tol_wall);
  //fs.exchange_all();
  //boil::plot->plot(fs,scp, "fs-clr", 0);
 
  return;
}

/***********************
 * ancillary function
 ***********************/
real Heaviside::fs_val(const Comp m, const int i, const int j, const int k,
                       const Scalar & nx, const Scalar & ny,
                       const Scalar & nz, const Scalar & nalpha) {

  real alpha = nalpha[i][j][k];
  /* degenerate case */
  if(!boil::realistic(alpha))
    return boil::unreal;

  real xpos = nalpha.xc(i);
  real ypos = nalpha.yc(j);
  real zpos = nalpha.zc(k);

  real vx = nx[i][j][k]; 
  real vy = ny[i][j][k]; 
  real vz = nz[i][j][k]; 

  if(m==Comp::i()) {
    if(fabs(vx)<boil::atto)
      return boil::unreal;
    else
      return (alpha-vy*ypos-vz*zpos)/vx;
  } else if(m==Comp::j()) {
    if(fabs(vy)<boil::atto)
      return boil::unreal;
    else
      return (alpha-vx*xpos-vz*zpos)/vy;
  } else {
    if(fabs(vz)<boil::atto)
      return boil::unreal;
    else
      return (alpha-vx*xpos-vy*ypos)/vz;
  }

  return boil::unreal;
}
