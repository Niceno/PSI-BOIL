#include "tif.h"

/***************************************************************************//**
 *  Returns the value of the interface temperature field
******************************************************************************/
real TIF::Tint(const int dir, const Comp &mcomp, real frac,
               const int i, const int j, const int k) const {
  if(variable_tif) {
    int of(1),off(0);
    if (dir < 0) {of = -1; off = -1;}

    if (mcomp == Comp::i()) {
      if((frac*tif.dxe(i+off)) < tif.dxc(i)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tif[i][j][k]+frac*tif[i+of][j][k];
    } else if (mcomp == Comp::j()) {
      if((frac*tif.dyn(j+off)) < tif.dyc(j)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tif[i][j][k]+frac*tif[i][j+of][k];
    } else { 
      if((frac*tif.dzt(k+off)) < tif.dzc(k)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tif[i][j][k]+frac*tif[i][j][k+of];
    }
  } else return tr;
}

real TIF::Tint_old(const int dir, const Comp &mcomp, real frac,
                   const int i, const int j, const int k) const {
  if(variable_tif) {
    int of(1),off(0);
    if (dir < 0) {of = -1; off = -1;}

    if (mcomp == Comp::i()) {
      if((frac*tif.dxe(i+off)) < tif.dxc(i)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tifold[i][j][k]+frac*tifold[i+of][j][k];
    } else if (mcomp == Comp::j()) {
      if((frac*tif.dyn(j+off)) < tif.dyc(j)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tifold[i][j][k]+frac*tifold[i][j+of][k];
    } else { 
      if((frac*tif.dzt(k+off)) < tif.dzc(k)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tifold[i][j][k]+frac*tifold[i][j][k+of];
    }
  } else return tr;
}

real TIF::Tint(const int i, const int j, const int k) const {
  if(variable_tif) 
    return tif[i][j][k];
  else
    return tr;
}

real TIF::Tint_old(const int i, const int j, const int k) const {
  if(variable_tif) 
    return tifold[i][j][k];
  else
    return tr;
}
