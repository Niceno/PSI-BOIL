#include "antoine.h"

/***************************************************************************//**
 *  Calculates interface temperature using the Antoine equation.
 *  IG approximation is used: alpha = p_alpha/p0.
 *
 *  Pressure at the interface assumed to be the same as the reference
 *  saturation pressure (corresponding to tsat).
 *
 *  Function output is stored in the tif field.
*******************************************************************************/
real Antoine::value(const int i, const int j, const int k) {
  /* epsilon is the vol. fraction of non-condensable gases */
  real alpha = 1.0-eps[i][j][k];

  if(alpha>0.0)
    return an_value(std::max(0.0,std::min(1.0,alpha)));
  else
    return 0.0;
}

real Antoine::an_value(const real alpha) {
  return B/(B/(C_K+tr)-log10(alpha))-C_K;
}

