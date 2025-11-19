#include "clapeyron.h"

/***************************************************************************//**
 *  Calculates interface temperature using the Clausius-Clapeyron IG
 *  approximation  with the latent heat taken from consistent (linear)
 *  formulation. Newton-Raphson method is used for iterations.
 *  cl and cv are assumed constant; tsati = 1/tsat.
 *
 *  Pressure at the interface assumed to be the same as the reference
 *  saturation pressure (corresponding to tsat).
 *
 *  errmax controls the N-R iterations on consistent solution.
 *
 *  Function output is stored in the tif field.
*******************************************************************************/
real Clapeyron::value(const int i, const int j, const int k) {
  /* epsilon is the vol. fraction of non-condensable gases */
  real alpha = 1.0-eps[i][j][k];

  if(alpha>0.0)
    return cc_iterate(std::max(0.0,std::min(1.0,alpha)));
  else
    return 0.0;
}

/* iteration */
real Clapeyron::cc_iterate(const real alpha) {
  /* initialisation [inconsistent c-c] */
  real tinv = tri - Rm/latent * log(alpha);

  real tau_m1 = tr*tinv;
  real tau_m0 = tau_m1;
  real t_err = err(tau_m0,alpha);

  int niter(0);
  while((fabs(t_err) > errmax) && niter != nmax) {
    niter++;
    tau_m0 = tau_m1 - err(tau_m1,alpha)/err_prime(tau_m1);
    t_err = err(tau_m0,alpha);
    tau_m1 = tau_m0;
  }

  return tr/tau_m0;
}

/* newton-raphson error evaluation functions, using temperature ratio tau
   and vapor volume fraction alpha */
real Clapeyron::err(const real tau, const real alpha) {
  return - log(alpha) 
         + latent_cst/(Rm*tr)*(1.0-tau)
         - latent_slp/Rm * log(tau);
}

real Clapeyron::err_prime(const real tau) {
  return -latent_cst/(Rm*tr) - latent_slp/Rm/tau;
}
