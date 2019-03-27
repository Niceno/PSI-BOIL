#include "vof.h"

/* liq+gas, both upw; note that jv is now positive always */
real VOF::iterate_flux(const real jv, const real dirsgn, const real cflrat,
                       const real alphaliq, const real alphagas,
                       const real vm1, const real vm2, const real vm3) {
  real x0 = 0.0;
  real x2 = std::min(1.0,jv/std::min(cflrat+boil::pico,1.0));

  real liq0 = calc_v_iter(alphaliq,vm1,vm2,vm3,x0,dirsgn);
  real liq2 = calc_v_iter(alphaliq,vm1,vm2,vm3,x2,dirsgn);
 
  real x0g = 0.0;
  real x2g = std::min(1.0,x2*cflrat);

  /* gas flows in the same direction and sense of n is inverted,
   * dirsgn -> -dirsgn */
  real gas0 = calc_v_iter(alphagas,vm1,vm2,vm3,x0g,-dirsgn);
  real gas2 = calc_v_iter(alphagas,vm1,vm2,vm3,x2g,-dirsgn);

  real j0 = liq0+gas0-jv;
  real j2 = liq2+gas2-jv;

  int iter(0);
  while(fabs(j0-j2)>tol_flux&&iter<maxiter) {
    iter++;
    
    /* midpoint */
    real x1 = 0.5*(x0+x2);
    real liq1 = calc_v_iter(alphaliq,vm1,vm2,vm3,x1,dirsgn);
    real x1g = std::min(1.0,x1*cflrat);
    real gas1 = calc_v_iter(alphagas,vm1,vm2,vm3,x1g,-dirsgn);
  
    real j1 = liq1+gas1-jv;

    /* ridders method */
    real sgnj0 = (j0>0.0) - (j0<0.0);

    real x3 = x1 + (x1-x0) * sgnj0*j1/ sqrt(j1*j1-j0*j2);

    real liq3 = calc_v_iter(alphaliq,vm1,vm2,vm3,x3,dirsgn);
    real x3g = std::min(1.0,x3*cflrat);
    real gas3 = calc_v_iter(alphagas,vm1,vm2,vm3,x3g,-dirsgn);

    real j3 = liq3+gas3-jv;

    if       (j3*j1<0.0) {
      x0 = x1;
      j0 = j1;
      x2 = x3;
      j2 = j3;
    } else if(j3*j0<0.0) {
      x2 = x3;
      j2 = j3;
    } else if(j3*j2<0.0) {
      x0 = x3;
      j0 = j3;
    } else if(j1*j2<0.0) {
      x0 = x1;
      j0 = j1;
    } else {
      x2 = x1;
      j2 = j1;
    }

    boil::oout<<"Iteration: "<<iter<<" | "<<x0<<" "<<x2<<" | "<<j0<<" "<<j2<<" | "<<x1<<" "<<x3<<" | "<<j1<<" "<<j3<<boil::endl;


  } /* iteration loop */

  /* interpolate bw j0 and j2 */
  real xtarget;
  if(fabs(j2-j0)<boil::pico) {
    xtarget = x0;
  } else {
    xtarget = -j0/(j2-j0) * (x2-x0) + x0;
  }

  return calc_v_iter(alphaliq,vm1,vm2,vm3,xtarget,dirsgn);
}

/* liq upw, gas downwind, jv can have either sign, cflrat is positive,
 * dirsgnup is positive if the normal vect upw pointing to liq aligns with jv
 * dirsgndn is negative if the normal vect dnw pointing to liq aligns with jv
 * alphaup is the liquid one, alphadn is the gas one,
 * liq contribution is positive, gas contribution is negative */
#if 1
real VOF::iterate_flux(const real jv, const real dirsgnup, const real dirsgndn,
                       const real cflrat, const real alphaup,
                       const real vm1up, const real vm2up, const real vm3up,
                       const real alphadn,
                       const real vm1dn, const real vm2dn, const real vm3dn,
                       const real x0start, const real x2start) {
  real x0 = x0start;
  real x2 = x2start;

  real liq0 = calc_v_iter(alphaup,vm1up,vm2up,vm3up,x0,dirsgnup);
  real liq2 = calc_v_iter(alphaup,vm1up,vm2up,vm3up,x2,dirsgnup);
 
  real x0g = std::min(1.0,x0*cflrat);
  real x2g = std::min(1.0,x2*cflrat);

  /* gas flows in the opposite direction and sense of n is inverted,
   * dirsgn -> -dirsgn and minus in front */
  real gas0 = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,x0g,-dirsgndn);
  real gas2 = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,x2g,-dirsgndn);

  real j0 = liq0+gas0-jv;
  real j2 = liq2+gas2-jv;

  assert(j0*j2<=0.0);

  int iter(0);
  while(fabs(j0-j2)>tol_flux&&iter<maxiter) {
    iter++;
    
    /* midpoint */
    real x1 = 0.5*(x0+x2);
    real liq1 = calc_v_iter(alphaup,vm1up,vm2up,vm3up,x1,dirsgnup);
    real x1g = std::min(1.0,x1*cflrat);
    real gas1 = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,x1g,-dirsgndn);
  
    real j1 = liq1+gas1-jv;

    /* ridders method */
    real sgnj0 = (j0>0.0) - (j0<0.0);

    real x3 = x1 + (x1-x0) * sgnj0*j1/ sqrt(j1*j1-j0*j2);

    real liq3 = calc_v_iter(alphaup,vm1up,vm2up,vm3up,x3,dirsgnup);
    real x3g = std::min(1.0,x3*cflrat);
    real gas3 = -calc_v_iter(alphadn,vm1dn,vm2dn,vm3dn,x3g,-dirsgndn);

    real j3 = liq3+gas3-jv;

    if       (j3*j1<0.0) {
      x0 = x1;
      j0 = j1;
      x2 = x3;
      j2 = j3;
    } else if(j3*j0<0.0) {
      x2 = x3;
      j2 = j3;
    } else if(j3*j2<0.0) {
      x0 = x3;
      j0 = j3;
    } else if(j1*j2<0.0) {
      x0 = x1;
      j0 = j1;
    } else {
      x2 = x1;
      j2 = j1;
    }

    //boil::oout<<"Iteration: "<<iter<<" | "<<x0<<" "<<x2<<" | "<<j0<<" "<<j2<<" | "<<x1<<" "<<x3<<" | "<<j1<<" "<<j3<<boil::endl;


  } /* iteration loop */

  /* interpolate bw j0 and j2 */
  real xtarget;
  if(fabs(j2-j0)<boil::pico) {
    xtarget = x0;
  } else {
    xtarget = -j0/(j2-j0) * (x2-x0) + x0;
  }

  return calc_v_iter(alphaup,vm1up,vm2up,vm3up,xtarget,dirsgnup);
}
#endif

real VOF::calc_v_iter(const real alpha, const real vm1, const real vm2,
                      const real vm3, const real absg, const real dirsgn) {
  real alp = alpha;
  real v1 = vm1;
  real ra = vm1 * (1.0 - absg);
  real qa = 1.0/(1.0-ra);
  if (dirsgn > 0) alp = alp -ra;
  v1 = v1 * absg;

  /* calculate f: reduced flux (for unit volume) */
  return calc_v(alp*qa, v1*qa, vm2*qa, vm3*qa) * absg;
}
