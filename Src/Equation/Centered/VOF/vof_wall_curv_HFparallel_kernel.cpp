#include "vof.h"

/******************************************************************************/
real VOF::wall_curv_HFparallel_kernel(const real hc, const real hp,
                                      const real dc, const real dp,
                                      const real mult,
                                      const real cang) {
/***************************************************************************//**
*  \brief Calculate curvature using hybrid height-function/divergence-of-normal
*         approach in 2D geometry, heights constructed in parallel direction.
*
*     output: kappa
*******************************************************************************/
#if 1
  real hzp = (hp-hc)/dp;

  real nzp = hzp/sqrt(1.+hzp*hzp)*-mult;
  real nxp = -1./sqrt(1.+hzp*hzp)*-mult;

  real nzm = cos(cang*mult);
  real nxm = sin(cang*mult);
  
  //boil::oout<<"wall: "<<hzp<<" "<<cang*180./boil::pi<<" "<<hc<<" "<<(dc+0.382331)/1.47721-nzp<<" | "<<nzp<<" "<<nxp<<" | "<<nzm<<" "<<nxm<<" | "<<(nzp-nzm)/dc/(1.0/1.72237)-1.<<" "<<0.5*(nxp+nxm)/hc/(1.0/1.72237)-1.<<boil::endl;

  return (nzp-nzm)/dc;
#else
  real dm = dc;
  real hm = hc + dm*mult*cos(cang)/sin(cang);
  real h_1c = (hp-hm)/(dm+dp);
  real h_1u = (hp-hc)/dp;
  real h_1d = (hc-hm)/dm;
  real h_11 = 2.*(dp*hm+dm*hp-(dp+dm)*hc)/dp/dm/(dp+dm);

  //boil::oout<<"wall: "<<hc<<" "<<hm<<" | "<<h_1u<<" "<<h_1c<<" "<<h_1d<<" | "<<h_11<<" "<<-mult*h_11/pow(1.+h_1c*h_1c,1.5)<<boil::endl;

  return -mult*h_11/pow(1.+h_1c*h_1c,1.5);
#endif
}

/******************************************************************************/
real VOF::wall_curv_HFparallel_kernel(const real hm,
                                      const real hc, const real hp,
                                      const real dm, const real dc,
                                      const real dp,
                                      const real mult) {
/***************************************************************************//**
*  \brief Calculate curvature using hybrid height-function/divergence-of-normal
*         approach in 2D geometry, heights constructed in parallel direction.
*
*     output: kappa
*******************************************************************************/

#if 0
  real hzm = (hc-hm)/dm;
  real hzp = (hp-hc)/dp;

  real nzp = hzp/sqrt(1.+hzp*hzp)*-mult;
  real nxp = -1./sqrt(1.+hzp*hzp)*-mult;

  real nzm = hzm/sqrt(1.+hzm*hzm)*-mult;
  real nxm = -1./sqrt(1.+hzm*hzm)*-mult;
  
  //boil::oout<<hzm<<" "<<hzp<<" | "<<nzp<<" "<<nxp<<" | "<<nzm<<" "<<nxm<<boil::endl;

  return (nzp-nzm)/dc;
#else 
  real h_1c = (hp-hm)/(dm+dp);
  real h_1u = (hp-hc)/dp;
  real h_1d = (hc-hm)/dm;
  real h_11 = 2.*(dp*hm+dm*hp-(dp+dm)*hc)/dp/dm/(dp+dm);

  //boil::oout<<"near: "<<hc<<" "<<hm<<" | "<<h_1u<<" "<<h_1c<<" "<<h_1d<<" | "<<h_11<<" "<<-mult*h_11/pow(1.+h_1c*h_1c,1.5)<<boil::endl;
  return -mult*h_11/pow(1.+h_1c*h_1c,1.5);
#endif
}
