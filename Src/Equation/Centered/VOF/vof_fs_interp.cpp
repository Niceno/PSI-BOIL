#include "vof.h"

/******************************************************************************/
void VOF::cal_fs_interp() {
/***************************************************************************//**
 \brief Calculate free-surface position between cell centers
    if there is no interface in the cell, unreal=yotta (=1e+24) is stored.
    output: fs
*******************************************************************************/

  /* initialize */
  for_m(m)
    for_avmijk(fs,m,i,j,k)
      fs[m][i][j][k] = boil::unreal;

  Comp m;

  /******************************************
  *             x-direction                 *
  ******************************************/

  m = Comp::i();
  for_wvmijk(fs,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */

    /* degenerate cases */
    real clrw = phi[i-1][j][k];
    real clre = phi[i  ][j][k];

    if((clrw-phisurf)*(clre-phisurf)>0.0)
      continue;

    if(  (clrw<boil::pico&&clre-1.0>-boil::pico)
       ||(clrw-1.0>-boil::pico&&clre<boil::pico)) {
      fs[m][i][j][k] = phi.xn(i);
      continue;
    }
 
    if(fabs(clrw-clre)<boil::pico) {
      fs[m][i][j][k] = phi.xn(i);
      continue;
    }

    fs[m][i][j][k] =  phi.xc(i-1) + (phisurf-clrw)/(clre-clrw)*phi.dxe(i-1);
  }

  /******************************************
  *             y-direction                 *
  ******************************************/

  m = Comp::j();
  for_wvmijk(fs,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */

    /* degenerate cases */
    real clrs = phi[i][j-1][k];
    real clrn = phi[i][j  ][k];

    if((clrs-phisurf)*(clrn-phisurf)>0.0)
      continue;

    if(  (clrs<boil::pico&&clrn-1.0>-boil::pico)
       ||(clrs-1.0>-boil::pico&&clrn<boil::pico)) {
      fs[m][i][j][k] = phi.yn(j);
      continue;
    }

    if(fabs(clrs-clrn)<boil::pico) {
      fs[m][i][j][k] = phi.yn(j);
      continue;
    }

    fs[m][i][j][k] =  phi.yc(j-1) + (phisurf-clrs)/(clrn-clrs)*phi.dyn(j-1);
  }

  /******************************************
  *             z-direction                 *
  ******************************************/

  m = Comp::k();
  for_wvmijk(fs,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */

    /* degenerate cases */
    real clrb = phi[i][j][k-1];
    real clrt = phi[i][j][k  ];

    if((clrb-phisurf)*(clrt-phisurf)>0.0)
      continue;

    if(  (clrb<boil::pico&&clrt-1.0>-boil::pico)
       ||(clrb-1.0>-boil::pico&&clrt<boil::pico)) {
      fs[m][i][j][k] = phi.zn(k);
      continue;
    }

    if(fabs(clrb-clrt)<boil::pico) {
      fs[m][i][j][k] = phi.zn(k);
      continue;
    }

    fs[m][i][j][k] =  phi.zc(k-1) + (phisurf-clrb)/(clrt-clrb)*phi.dzt(k-1);
  }

  /* correct at boundaries */
  //fs_bnd();
  //fs.exchange_all();

  //boil::plot->plot(fs,phi, "fs-clr", 0);

  return;
}
