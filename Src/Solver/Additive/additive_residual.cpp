#include "additive.h"

/******************************************************************************/
real AC::residual(Centered & h, real * linf) const {
	
  /* estimate residual */
  h.res = h.fnew - h.A * h.phi;
  real r2 = h.res.dot(h.res);

  /* this will be used for normalisation */
  real f2 = h.fnew.dot(h.fnew);

  if(linf && r2>0.)
    *linf = h.res.max_abs()/sqrt(r2);

  return sqrt(r2 / (f2+boil::pico));
}	
