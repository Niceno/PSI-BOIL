#include "vof.h"

/******************************************************************************/
void VOF::cal_fs3() {
/***************************************************************************//**
 \brief Calculate free-surface position between cell centers
    if there is no interface in the cell, yotta (=1e+24) is stored.
    plane: vm1*x + vm2*y + vm3*z = alpha
    output: fs
*******************************************************************************/

  /* tolerance is necessary because the linear-planes are not closed */
  //real tol = 0.0e-3; 
  //real tol = 2.0e-2;  
  real tol = 0.0e-2;
        
  /* initialize */
  for_m(m)
    for_avmijk(fs,m,i,j,k)
      fs[m][i][j][k] = boil::yotta;

  Comp m;
  
  /******************************************
  *             x-direction                 *
  ******************************************/

  m = Comp::i();
  //for_vmijk(u,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */
  for(int i=si(); i<=ei()+1; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()  ; k++) {
  
#if 0
    if(i<6&&j==1&&k==1)
      boil::oout<<"VOF::FS3 "<<i<<" "<<iflag[i-1][j][k]<<" "<<iflag[i][j][k]<<boil::endl;  
#endif
             
    //if(abs(iflag[i-1][j][k])>=1||abs(iflag[i][j][k])>=1)
    //  continue;

    /* degenerate cases */
    real clrw = phi[i-1][j][k];
    real clre = phi[i  ][j][k];

    if((clrw-phisurf)*(clre-phisurf)>0.0) 
      continue;
#if 0
    if(i<6&&j==1&&k==1)
      boil::oout<<"           "<<clrw<<" "<<clre<<boil::endl;  
#endif  
  
    if(  (clrw<boil::pico&&clre>1.0-boil::pico)
       ||(clrw>1.0-boil::pico&&clre<boil::pico)) {
      fs[m][i][j][k] = phi.xn(i);
      continue;
    }

    /* calculate normalized candidate positions */
    real fsxw = fs_val(m,i-1,j,k);
    real fsxe = fs_val(m,i  ,j,k);

#if 0
    if(i<6&&j==1&&k==1)
      boil::oout<<"           "<<fsxw<<" "<<fsxe<<boil::endl;  
#endif

    /* --------------------------------------------
       first check: are candidate positions 
       a. between cell centres
       b. inside the correct cell (real interface)
    ----------------------------------------------- */
    bool flagw = (0.5     <= fsxw && fsxw <= 1.0+tol);
    bool flage = (0.0-tol <= fsxe && fsxe <= 0.5    );

    if     ( flagw && !flage) { /* west is real and east is not */
      fs[m][i][j][k] = phi.xn(i-1) + phi.dxc(i-1) * fsxw;

#if 0
      if(i<6&&j==1&&k==1)
        boil::oout<<"           "<<flagw<<" "<<flage<<" "<<fs[m][i][j][k]<<boil::endl;  
#endif
      continue;
    }
    else if(!flagw &&  flage) { /* east is real and west is not */
      fs[m][i][j][k] = phi.xn(i) + phi.dxc(i) * fsxe;

#if 0
      if(i<6&&j==1&&k==1)
        boil::oout<<"           "<<flagw<<" "<<flage<<" "<<fs[m][i][j][k]<<boil::endl;  
#endif
      continue;
    }
    else if( flagw &&  flage) { /* both are real */
#if 0
      fs[m][i][j][k] = phi.xn(i);
#else
      /* calculate distance between the two candidate positions */
      real fsx_diff = fsxe + 1.0 - fsxw;
      if(fabs(fsx_diff)<tol) { /* they are close to each other */
        fs[m][i][j][k] = phi.xn(i);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = (1.0-fsxw>fsxe) ? 
                         phi.xn(i-1) + phi.dxc(i-1) * fsxw
                       : phi.xn(i  ) + phi.dxc(i  ) * fsxe;
      }
#endif

#if 0
      if(i<6&&j==1&&k==1)
        boil::oout<<"           "<<flagw<<" "<<flage<<" "<<fs[m][i][j][k]<<boil::endl;  
#endif
      continue;
    }
#if 0
    /* -----------------------------------------------
       second check: are candidate positions 
       a. between cell centres
       b. inside the wrong cell (imaginary interface)
    -------------------------------------------------- */
    flagw = ( 1.0-tol <= fsxw && fsxw <= 1.5    );
    flage = (-0.5     <= fsxe && fsxe <= 0.0+tol);
    if     ( flagw && !flage) { /* west is imaginary and east is not */
      fs[m][i][j][k] = phi.xn(i-1) + phi.dxc(i-1) * fsxw;
      continue;
    }
    else if(!flagw &&  flage) { /* east is imaginary and west is not */
      fs[m][i][j][k] = phi.xn(i) + phi.dxc(i) * fsxe;
      continue;
    }
    else { /* both are either imaginary or not intersecting */
      fs[m][i][j][k] = phi.xn(i);
      continue;
    }
#else
    fs[m][i][j][k] = phi.xn(i);
#endif
  } 
  
  /******************************************
  *             y-direction                 *
  ******************************************/

  m = Comp::j();
  //for_vmijk(u,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()+1; j++)
  for(int k=sk(); k<=ek()  ; k++) {
               
    //if(abs(iflag[i][j-1][k])>=1||abs(iflag[i][j][k])>=1)
    //  continue;

    /* degenerate cases */
    real clrs = phi[i][j-1][k];
    real clrn = phi[i][j  ][k];

    if((clrs-phisurf)*(clrn-phisurf)>0.0) 
      continue;
    
    if(  (clrs<boil::pico&&clrn>1.0-boil::pico)
       ||(clrs>1.0-boil::pico&&clrn<boil::pico)) {
      fs[m][i][j][k] = phi.yn(j);
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
    bool flags = (0.5     <= fsys && fsys <= 1.0+tol);
    bool flagn = (0.0-tol <= fsyn && fsyn <= 0.5+tol);
   
    if     ( flags && !flagn) { /* south is real and north is not */
      fs[m][i][j][k] = phi.yn(j-1) + phi.dyc(j-1) * fsys;

     //if(time->current_step()==518&&i==5)
     //  boil::oout<<i<<" "<<j<<" "<<k<<" "<<fs[m][i][j][k]<<" "<<flags<<" "<<flagn<<" "<<fsys<<" "<<fsyn<<" "<<phi.yn(j-1) + phi.dyc(j-1) * fsys<<" "<<phi.yn(j) + phi.dyc(j) * fsyn<<boil::endl;
      continue;
    }
    else if(!flags &&  flagn) { /* north is real and south is not */
      fs[m][i][j][k] = phi.yn(j) + phi.dyc(j) * fsyn;

     //if(time->current_step()==518&&i==5)
     //   boil::oout<<i<<" "<<j<<" "<<k<<" "<<fs[m][i][j][k]<<" "<<flags<<" "<<flagn<<" "<<fsys<<" "<<fsyn<<" "<<phi.yn(j-1) + phi.dyc(j-1) * fsys<<" "<<phi.yn(j) + phi.dyc(j) * fsyn<<boil::endl;
      continue;
    }
    else if( flags &&  flagn) { /* both are real */
#if 0
      fs[m][i][j][k] = phi.yn(j);
#else
      /* calculate distance between the two candidate positions */
      real fsy_diff = fsyn + 1.0 - fsys;
      if(fabs(fsy_diff)<tol) { /* they are close to each other */
        fs[m][i][j][k] = phi.yn(j);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = (1.0-fsys>fsyn) ? 
                         phi.yn(j-1) + phi.dyc(j-1) * fsys
                       : phi.yn(j  ) + phi.dyc(j  ) * fsyn;
      }

     //if(time->current_step()==518&&i==5)
     //   boil::oout<<i<<" "<<j<<" "<<k<<" "<<fs[m][i][j][k]<<" "<<flags<<" "<<flagn<<" "<<fsys<<" "<<fsyn<<" "<<phi.yn(j-1) + phi.dyc(j-1) * fsys<<" "<<phi.yn(j) + phi.dyc(j) * fsyn<<boil::endl;
#endif
    
      continue;
    }
#if 0     
    /* -----------------------------------------------
       second check: are candidate positions 
       a. between cell centres
       b. inside the wrong cell (imaginary interface)
    -------------------------------------------------- */
    flags = ( 1.0-tol <= fsys && fsys <= 1.5    );
    flagn = (-0.5     <= fsyn && fsyn <= 0.0+tol);
    if     ( flags && !flagn) { /* south is imaginary and north is not */
      fs[m][i][j][k] = phi.yn(j-1) + phi.dyc(j-1) * fsys;
      continue;
    }
    else if(!flags &&  flagn) { /* north is imaginary and south is not */
      fs[m][i][j][k] = phi.yn(j) + phi.dyc(j) * fsyn;
      continue;
    }
    else { /* both are either imaginary or not intersecting */
      fs[m][i][j][k] = phi.yn(j);
      continue;
    }
#else
    fs[m][i][j][k] = phi.yn(j);

    //if(time->current_step()==518&&i==5)
    //    boil::oout<<i<<" "<<j<<" "<<k<<" "<<fs[m][i][j][k]<<" "<<flags<<" "<<flagn<<" "<<fsys<<" "<<fsyn<<" "<<phi.yn(j-1) + phi.dyc(j-1) * fsys<<" "<<phi.yn(j) + phi.dyc(j) * fsyn<<boil::endl;
#endif
  } 
  
  /******************************************
  *             z-direction                 *
  ******************************************/

  m = Comp::k();
  //for_vmijk(u,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()+1; k++) {
               
    //if(abs(iflag[i][j][k-1])>=1||abs(iflag[i][j][k])>=1)
    //  continue;

    /* degenerate cases */
    real clrb = phi[i][j][k-1];
    real clrt = phi[i][j][k  ];

    if((clrb-phisurf)*(clrt-phisurf)>0.0) 
      continue;
    
    if(  (clrb<boil::pico&&clrt>1.0-boil::pico)
       ||(clrb>1.0-boil::pico&&clrt<boil::pico)) {
      fs[m][i][j][k] = phi.zn(k);
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
    bool flagb = (0.5     <= fszb && fszb <= 1.0+tol);
    bool flagt = (0.0-tol <= fszt && fszt <= 0.5    );
   
    if     ( flagb && !flagt) { /* south is real and north is not */
      fs[m][i][j][k] = phi.zn(k-1) + phi.dzc(k-1) * fszb;
      continue;
    }
    else if(!flagb &&  flagt) { /* north is real and south is not */
      fs[m][i][j][k] = phi.zn(k) + phi.dzc(k) * fszt;
      continue;
    }
    else if( flagb &&  flagt) { /* both are real */
#if 0
      fs[m][i][j][k] = phi.zn(k);
#else
      /* calculate distance between the two candidate positions */
      real fsz_diff = fszt + 1.0 - fszb;
      if(fabs(fsz_diff)<tol) { /* they are close to each other */
        fs[m][i][j][k] = phi.zn(k);
      } else {  /* choose the "more real" value */
        fs[m][i][j][k] = (1.0-fszb>fszt) ? 
                         phi.zn(k-1) + phi.dzc(k-1) * fszb
                       : phi.zn(k  ) + phi.dzc(k  ) * fszt;
      }
#endif
      continue;
    }
#if 0     
    /* -----------------------------------------------
       second check: are candidate positions 
       a. between cell centres
       b. inside the wrong cell (imaginary interface)
    -------------------------------------------------- */
    flagb = ( 1.0-tol <= fszb && fszb <= 1.5    );
    flagt = (-0.5     <= fszt && fszt <= 0.0+tol);
    if     ( flagb && !flagt) { /* south is imaginary and north is not */
      fs[m][i][j][k] = phi.zn(k-1) + phi.dzc(k-1) * fszb;
      continue;
    }
    else if(!flagb &&  flagt) { /* north is imaginary and south is not */
      fs[m][i][j][k] = phi.zn(k) + phi.dzc(k) * fszt;
      continue;
    }
    else { /* both are either imaginary or not intersecting */
      fs[m][i][j][k] = phi.zn(k);
      continue;
    }
#else
    fs[m][i][j][k] = phi.zn(k);
#endif
  } 

  //fs.exchange_all();

  //boil::plot->plot(fs,phi, "fs-clr", 0);
 
  return;
}

/***********************
 * ancillary function
 ***********************/
real VOF::fs_val(const Comp m, const int i, const int j, const int k) {

  /* calculate vn1, vn2, vn3: normal vector at face center */
  real vn1 = -nx[i][j][k];
  real vn2 = -ny[i][j][k];
  real vn3 = -nz[i][j][k];

  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real qa = 1.0/(vm1+vm2+vm3);
  vm1 *= qa;
  vm2 *= qa;
  vm3 *= qa;
  real c = phi[i][j][k];
  real alpha = calc_alpha(c, vm1, vm2, vm3);
  
#if 0
  if(i<6&&j==1&&k==1)
    boil::oout<<"VOF::FSV "<<i<<" "<<vn1<<" "<<vn2<<" "<<vn3<<" "<<alpha<<boil::endl;  
#endif

  real xpos = 0.5;
  real ypos = 0.5;
  real zpos = 0.5;

  if(m==Comp::i()) {
    real xuni = (alpha-vm2*ypos-vm3*zpos)/(vm1+boil::pico);
    if(vn1<0)
      xuni = 1.0-xuni;
    return xuni;
    //return phi.xn(i) + phi.dxc(i) * xuni;
  } else if(m==Comp::j()) {
    real yuni = (alpha-vm1*xpos-vm3*zpos)/(vm2+boil::pico);
    if(vn2<0)
      yuni = 1.0-yuni;
    return yuni;
    //return phi.yn(j) + phi.dyc(j) * yuni;
  } else {
    real zuni = (alpha-vm1*xpos-vm2*ypos)/(vm3+boil::pico);
    if(vn3<0)
      zuni = 1.0-zuni;
    return zuni;
    //return phi.zn(k) + phi.dzc(k) * zuni;
  }

  return 0.0;
}

