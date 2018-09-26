#include "scalar.h"

/******************************************************************************/
real Scalar::min() const {
  real vmin = val[1][1][1];
  for_ijk(i,j,k) {
    if( val[i][j][k] < vmin ) {
      vmin = val[i][j][k];
    }
  }
  boil::cart.min_real(&vmin);
  return vmin;
}

/******************************************************************************/
real Scalar::min_abs() const {
  real vmin = fabs(val[1][1][1]);
  for_ijk(i,j,k) {
    if( fabs(val[i][j][k]) < vmin ) {
      vmin = fabs(val[i][j][k]);
    }
  }
  boil::cart.min_real(&vmin);
  return vmin;
}

/******************************************************************************/
real Scalar::min_at() const {
  real vmin = val[1][1][1];
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k] < vmin ) {
      vmin = val[i][j][k]; im=i; jm=j; km=k;
    }
  }
  real vmin_l = vmin;
  boil::cart.min_real(&vmin);
  if(vmin == vmin_l) {
    boil::aout << name() << " min = " << vmin << " at: " 
               << xc(im) << " " 
               << yc(jm) << " " 
               << zc(km) << boil::endl;
  }
  return vmin;
}

/******************************************************************************/
real Scalar::max() const {
  real vmax = val[1][1][1];
  for_ijk(i,j,k) {
    if( val[i][j][k] > vmax ) {
      vmax = val[i][j][k];
    }
  }
  boil::cart.max_real(&vmax);
  return vmax;
}

/******************************************************************************/
real Scalar::max_abs() const {
  real vmax = fabs(val[1][1][1]);
  for_ijk(i,j,k) {
    if( fabs(val[i][j][k]) > vmax ) {
      vmax = fabs(val[i][j][k]);
    }
  }
  boil::cart.max_real(&vmax);
  return vmax;
}

/******************************************************************************/
real Scalar::max_at() const {
  real vmax = val[1][1][1];
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k] > vmax ) {
      vmax = val[i][j][k]; im=i; jm=j; km=k;
    }
  }
  real vmax_l = vmax;
  boil::cart.max_real(&vmax);
  if(vmax == vmax_l) {
    boil::aout << name() << " max = " << vmax << " at: " 
               << xc(im) << " " 
               << yc(jm) << " " 
               << zc(km) << boil::endl;
  }
  return vmax;
}

/*-----------------------------------------------------------------------------+
 '$Id: scalar_minmax.cpp,v 1.7 2011/09/23 12:32:17 niceno Exp $'/
+-----------------------------------------------------------------------------*/
