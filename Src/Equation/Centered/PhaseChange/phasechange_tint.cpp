#include "phasechange.h"

/***************************************************************************//**
 *  Returns interface temperature (calculated in EnthalpyTIF) 
*******************************************************************************/
real PhaseChange::Tint(const int i, const int j, const int k) {
  if (tif) return (*tif)[i][j][k];
  else return tsat;
}

real PhaseChange::Tint(const int dir, const Comp &mcomp, real frac,
                      const int i, const int j, const int k) {
  if (tif) {
    int of(1),off(0);
    if (dir < 0) {of = -1; off = -1;}
 
    if (mcomp == Comp::i()) {
      if((frac*phi.dxe(i+off)) < phi.dxc(i)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*(*tif)[i][j][k]+frac*(*tif)[i+of][j][k];
    } else if (mcomp == Comp::j()) {
      if((frac*phi.dyn(j+off)) < phi.dyc(j)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*(*tif)[i][j][k]+frac*(*tif)[i][j+of][k];
    } else {
      if((frac*phi.dzt(k+off)) < phi.dzc(k)/2.0) frac = 0.0; else frac = 1.0;
      return (1.0-frac)*(*tif)[i][j][k]+frac*(*tif)[i][j][k+of];
    }
  } else return tsat;
}

