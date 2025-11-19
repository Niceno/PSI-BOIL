#include "antoine.h" 

/* forward function */
real Antoine::temperature(const real eps) {
  if(eps<1.0)
    return an_value(std::max(0.0,std::min(1.0,1.0-eps)));
  else
    return 0.0;
}

/* inverse function */
real Antoine::epsilon(const real tint) {
  if((tint+C_K)<=0.)
    return 1.0;
  else
    return 1.-pow(10.,B/(C_K+tr)-B/(C_K+tint));
}

/* saturation pressure */
real Antoine::pressure(const real tint) {
  /* log10pressure in mmHg */
  real log10phg = A - B/(C_K+tint);

  /* to Pa */
  return 133.322*pow(10.,log10phg);
}
