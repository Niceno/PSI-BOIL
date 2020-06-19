#include "vof.h"

/******************************************************************************/
real VOF::wall_curv_HFnormal_kernel(const real x0,
                                    const real hm, const real hc, const real hp,
                                    const real dm, const real dc, const real dp,
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
    //nxm = -mult*sin(cang);
    nxm = -sin(cang);
    nzm =  cos(cang);
  }

  if(hp>0.) {
    real hxp = (hp-hc)/dp;
    nxp = hxp/sqrt(1.+hxp*hxp)*-mult;
    nzp = -1./sqrt(1.+hxp*hxp)*-mult;
  } else {
    //nxp =  mult*sin(cang);
    nxp =  sin(cang);
    nzp =  cos(cang);
  }

  //boil::oout<<nxp<<" "<<nxm<<" "<<dc<<boil::endl;

  /* droplets have positive curvature */
  return (nxp-nxm)/dc;
#else
  real hhm = hm;
  real hhp = hp;
  if(hm<=0.) {
    hhm = hc - dm*tan(mult*cang);
  }
  if(hp<=0.) {
    hhp = hc - dc*tan(mult*cang);
  }

  real h_1c = (hhp-hhm)/(dm+dp);
  real h_1u = (hhp-hc)/dp;
  real h_1d = (hc-hhm)/dm;
  real h_11 = 2.*(dp*hhm+dm*hhp-(dp+dm)*hc)/dp/dm/(dp+dm);

  boil::oout<<hhm<<" "<<hc<<" "<<hhp<<" "<<dc<<boil::endl;

  return -mult*h_11/pow(1.+h_1c*h_1c,1.5);
#endif
}
