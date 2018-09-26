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

/*-----------------------------------------------------------------------------+
 '$Id: global_approx.cpp,v 1.4 2014/02/03 11:09:19 sato Exp $'/
+-----------------------------------------------------------------------------*/
