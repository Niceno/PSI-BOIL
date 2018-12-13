#include "enthalpytif.h"

/***************************************************************************//**
 *  Call to tifmodel
******************************************************************************/
real EnthalpyTIF::Tint(const int dir, const Comp &mcomp, real frac,
                      const int i, const int j, const int k) {
  return tifmodel.Tint(dir,mcomp,frac,i,j,k);
}

real EnthalpyTIF::Tint_old(const int dir, const Comp &mcomp, real frac,
                          const int i, const int j, const int k) {
  return tifmodel.Tint_old(dir,mcomp,frac,i,j,k);
}

real EnthalpyTIF::Tint(const int i, const int j, const int k) {
  return tifmodel.Tint(i,j,k);
}

real EnthalpyTIF::Tint_old(const int i, const int j, const int k) {
  return tifmodel.Tint_old(i,j,k);
}
