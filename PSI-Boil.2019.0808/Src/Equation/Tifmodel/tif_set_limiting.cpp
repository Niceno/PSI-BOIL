#include "tif.h"

/***************************************************************************//**
 * Set a limiting method for tif
******************************************************************************/
void TIF::set_weak_limiting(const real Tmin, const real Tmax) {
  weaklim = true;
  tmin = Tmin;
  tmax = Tmax;

  return;
}
