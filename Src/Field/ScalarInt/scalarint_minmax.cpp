#include "scalarint.h"

/******************************************************************************/
int ScalarInt::min() const {
  int vmin = val[1][1][1];
  for_ijk(i,j,k) {
    if( val[i][j][k] < vmin ) {
      vmin = val[i][j][k];
    }
  }
  boil::cart.min_int(&vmin);
  return vmin;
}

/******************************************************************************/
int ScalarInt::min_abs() const {
  int vmin = abs(val[1][1][1]);
  for_ijk(i,j,k) {
    if( abs(val[i][j][k]) < vmin ) {
      vmin = abs(val[i][j][k]);
    }
  }
  boil::cart.min_int(&vmin);
  return vmin;
}

/******************************************************************************/
int ScalarInt::min_at() const {
  int vmin = val[1][1][1];
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k] < vmin ) {
      vmin = val[i][j][k]; im=i; jm=j; km=k;
    }
  }
  int vmin_l = vmin;
  boil::cart.min_int(&vmin);
  if(vmin == vmin_l) {
    boil::aout << name() << " min = " << vmin << " at: " 
               << xc(im) << " " 
               << yc(jm) << " " 
               << zc(km) << boil::endl;
  }
  return vmin;
}

/******************************************************************************/
int ScalarInt::max() const {
  int vmax = val[1][1][1];
  for_ijk(i,j,k) {
    if( val[i][j][k] > vmax ) {
      vmax = val[i][j][k];
    }
  }
  boil::cart.max_int(&vmax);
  return vmax;
}

/******************************************************************************/
int ScalarInt::max_abs() const {
  int vmax = abs(val[1][1][1]);
  for_ijk(i,j,k) {
    if( abs(val[i][j][k]) > vmax ) {
      vmax = abs(val[i][j][k]);
    }
  }
  boil::cart.max_int(&vmax);
  return vmax;
}

/******************************************************************************/
int ScalarInt::max_at() const {
  int vmax = val[1][1][1];
  int  im   = 0;
  int  jm   = 0;
  int  km   = 0;
  for_ijk(i,j,k) {
    if( val[i][j][k] > vmax ) {
      vmax = val[i][j][k]; im=i; jm=j; km=k;
    }
  }
  int vmax_l = vmax;
  boil::cart.max_int(&vmax);
  if(vmax == vmax_l) {
    boil::aout << name() << " max = " << vmax << " at: " 
               << xc(im) << " " 
               << yc(jm) << " " 
               << zc(km) << boil::endl;
  }
  return vmax;
}
