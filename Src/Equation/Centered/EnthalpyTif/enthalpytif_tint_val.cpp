#include "enthalpytif.h"

/***************************************************************************//**
 *  Returns the value of the interface temperature field
******************************************************************************/
real EnthalpyTIF::Tint(const int dir, const Comp &mcomp, real frac,
                      const int i, const int j, const int k) {
  if(pres||mflx) {
    int of(1),off(0);
    if (dir < 0) {of = -1; off = -1;}

    if (mcomp == Comp::i()) {
      if((frac*phi.dxe(i+off)) < phi.dxc(i)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tif[i][j][k]+frac*tif[i+of][j][k];
    } else if (mcomp == Comp::j()) {
      if((frac*phi.dyn(j+off)) < phi.dyc(j)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tif[i][j][k]+frac*tif[i][j+of][k];
    } else { 
      if((frac*phi.dzt(k+off)) < phi.dzc(k)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tif[i][j][k]+frac*tif[i][j][k+of];
    }
  } else return tsat;
}

real EnthalpyTIF::Tint_old(const int dir, const Comp &mcomp, real frac,
                          const int i, const int j, const int k) {
  if(pres||mflx) {
    int of(1),off(0);
    if (dir < 0) {of = -1; off = -1;}

    if (mcomp == Comp::i()) {
      if((frac*phi.dxe(i+off)) < phi.dxc(i)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tifold[i][j][k]+frac*tifold[i+of][j][k];
    } else if (mcomp == Comp::j()) {
      if((frac*phi.dyn(j+off)) < phi.dyc(j)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tifold[i][j][k]+frac*tifold[i][j+of][k];
    } else { 
      if((frac*phi.dzt(k+off)) < phi.dzc(k)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*tifold[i][j][k]+frac*tifold[i][j][k+of];
    }
  } else return tsat;
}

real EnthalpyTIF::Tint(const int i, const int j, const int k) {
  if(pres||mflx) 
    return tif[i][j][k];
  else
    return tsat;
}

real EnthalpyTIF::Tint_old(const int i, const int j, const int k) {
  if(pres||mflx) 
    return tifold[i][j][k];
  else
    return tsat;
}
