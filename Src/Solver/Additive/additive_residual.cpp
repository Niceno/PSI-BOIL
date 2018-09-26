#include "additive.h"

/******************************************************************************/
real AC::residual(Centered & h) const {
	
  /* estimate residual */
  h.res = h.fnew - h.A * h.phi;
  real r2 = h.res.dot(h.res);

  /* this will be used for normalisation */
  real f2 = h.fnew.dot(h.fnew);

  return sqrt(r2 / (f2+boil::pico));
}	

/*-----------------------------------------------------------------------------+
 '$Id: additive_residual.cpp,v 1.5 2011/09/30 14:27:30 niceno Exp $'/
+-----------------------------------------------------------------------------*/
