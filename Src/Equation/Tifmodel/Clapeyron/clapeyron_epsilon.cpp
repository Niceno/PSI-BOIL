#include "clapeyron.h" 

/* forward function */
real Clapeyron::temperature(const real eps) {
  if(eps<1.0)
    return cc_iterate(std::max(0.0,std::min(1.0,1.0-eps)));
  else
    return 0.0;
}

/* inverse function */
real Clapeyron::epsilon(const real tint) {
  if (tint > 0.0) {
    real logalp = latent_cst/(Rm*tr) * (1.0 - tr/tint)
                + latent_slp/Rm * log(tr/tint);
    return 1.0-exp(logalp);
  } else {
    return 1.0;
  }
}
