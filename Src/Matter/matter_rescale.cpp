#include "matter.h"
  
/******************************************************************************/
void Matter::rescale(const real xmult, const real tmult, const real mmult) {
/***************************************************************************//**
 \brief Rescale properties by length, time & mass coefficients.
        E.g. xmult = 1000 -> [m] -> [mm]
             tmult = 1e6  -> [s] -> [us]
             mmult = 1e-3 -> [kg]-> [ t]

  rho:   M/L^3 = kg/mmult/(m/xmult)^3 = kg/m^3 * mmult * xmult^3 
                      -> rho'   = rho    *mmult/xmult^3
  mu :   M/L/T        -> mu'    = mu     *mmult/xmult/tmult 
  cp :   M/T^2/L /K   -> cp'    = cp     *mmult/xmult/tmult^2
  lambda:M*L/T^3 /K   -> lambda'= lambda *mmult*xmult/tmult^3
  gamma: M/L/T =      -> gamma' = gamma  *mmult/xmult/tmult
  beta:  1/K          -> beta'  = beta   *1
  sigma: M/T^2 =      -> sigma' = sigma  *mmult/tmult/tmult

  furthermore, regarding observables

  kappa'    = kappa    *1/xmult              [1/L]
  press'    = press    *mmult/xmult/tmult^2  [M/L/T^2]
  velocity' = velocity *xmult/tmult          [L/T]
  force'    = force    *mmult*xmult/tmult^2  [M*L/T^2]
  masssrc'  = masssrc  *mmult/xmult^3/tmult  [M/L^3/T]

  choosing mmult = xmult*tmult and tmult = xmult results in

  mu'     = mu
  cp'     = cp/tmult
  lambda' = lambda
  gamma'  = gamma
  sigma'  = sigma

  as well as

  kappa'    = kappa    *1/xmult
  press'    = press    *1/tmult
  velocity' = velocity 
  force'    = force    *xmult
  masssrc'  = masssrc  *1/xmult^2

*******************************************************************************/

  boil::oout<<"Material properties adjusted by L,T,M factors "
            <<xmult<<" "<<tmult<<" "<<mmult
            <<boil::endl;

  if(!mixt) {
    rho(rho()->value()      *mmult/xmult/xmult/xmult);
    mu(mu()->value()        *mmult/xmult/tmult);
    cp(cp()->value()        *mmult/xmult/tmult/tmult);
    lambda(lambda()->value()*mmult*xmult/tmult/tmult/tmult);
    gamma(gamma()->value()  *mmult/xmult/tmult);
    beta(beta()->value());
  } else {
    sigma(sigma()->value()  *mmult/tmult/tmult);
  }

  return;
 
}
