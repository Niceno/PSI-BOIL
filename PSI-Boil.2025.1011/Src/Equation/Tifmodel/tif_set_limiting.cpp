#include "tif.h"

/***************************************************************************//**
 * Set a limiting method for tif
******************************************************************************/
void TIF::set_weak_limiting(const real Tmin, const real Tmax) {
  weaklim = true;
  //stronglim = false;
  tmin = Tmin;
  tmax = Tmax;

  return;
}

void TIF::set_strong_limiting(const Scalar * Tpr,
                              const Scalar * Clr,
                              const real Clrsurf) {
  //weaklim = false;
  stronglim = true;

  clrsurf = Clrsurf;
  clr = Clr;
  tpr = Tpr;

  return;
}

