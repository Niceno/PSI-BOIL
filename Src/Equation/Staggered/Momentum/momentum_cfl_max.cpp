#include "momentum.h"

/******************************************************************************/
real Momentum::cfl_max() const {
 
  real cfl, cflm = -1.0;

  Comp dir;
  int im  =  0;
  int jm  =  0;
  int km  =  0;

  Comp m = Comp::u();
  for_mijk(m,i,j,k) {
    cfl = fabs(u[m][i][j][k]) * time->dt() / u.dxc(m,i);
    if( cfl > cflm ) {cflm = cfl; dir = m; im=i; jm=j; km=k;} 
  }
  
  m = Comp::v();
  for_mijk(m,i,j,k) {
    cfl = fabs(u[m][i][j][k]) * time->dt() / u.dyc(m,j);
    if( cfl > cflm ) {cflm = cfl; dir = m; im=i; jm=j; km=k;}
  }

  m = Comp::w();
  for_mijk(m,i,j,k) {
    cfl = fabs(u[m][i][j][k]) * time->dt() / u.dzc(m,k);
    if( cfl > cflm ) {cflm = cfl; dir = m; im=i; jm=j; km=k;}
  }

  real cflm_l = cflm;
  boil::cart.max_real(&cflm);

  if(cflm == 0.0){
    boil::oout << "cfl max = " << cflm <<  boil::endl;
  } else if( approx(cflm,cflm_l) ){
    boil::aout << "@cfl_max; cfl max = " 
               << cflm 
               << " in direction " << dir 
               << " at: " 
               << u.xc(dir,im) << " " << u.yc(dir,jm) << " " << u.zc(dir,km) 
               << boil::endl;
  }

  return cflm;
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_cfl_max.cpp,v 1.15 2015/07/15 13:40:48 sato Exp $'/
+-----------------------------------------------------------------------------*/
