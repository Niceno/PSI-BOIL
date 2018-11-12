#include "global_precision.h"
#include "global_approx.h"

/******************************************************************************/
bool approx(real a, real b) {
  return approx(a, b, boil::pico);
}

/******************************************************************************/
bool approx(real a, real b, real diff) {

  if( fabs(a-b) < diff ) return true;
  else                   return false;
}
