#include "vofaxisym.h"

/******************************************************************************/
real VOFaxisym::wall_curv_HFnormal_kernel(const real x0, const real hm, 
                                          const real hc, const real hp,
                                          const real dm, 
                                          const real dc, const real dp,
                                          const real mult, const real cang) {
/***************************************************************************//**
*  \brief Calculate curvature using hybrid height-function/divergence-of-normal
*         approach in 2D geometry, heights constructed in normal direction.
*
*     output: kappa
*******************************************************************************/
#if 1
  /* we calculate normal vector pointing to the gas,
     heights are either gas or liquid, based on the mult sign */
  real nxm, nxp, nzm, nzp;
  
  if(hm>0.) {
    real hxm = (hc-hm)/dm;
    nxm = hxm/sqrt(1.+hxm*hxm)*-mult;
    nzm = -1./sqrt(1.+hxm*hxm)*-mult;
  } else {
    nxm = -sin(cang);
    nzm =  cos(cang);
  }

  if(hp>0.) {
    real hxp = (hp-hc)/dp;
    nxp = hxp/sqrt(1.+hxp*hxp)*-mult;
    nzp = -1./sqrt(1.+hxp*hxp)*-mult;
  } else {
    nxm =  sin(cang);
    nzm =  cos(cang);
  }

  /* droplets have positive curvature */
  return (nxp-nxm)/dc + 0.5*(nxp+nxm)/x0;
#else
  /* needs changing to normal */
  exit(0);
#endif
}
