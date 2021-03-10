#include "vof.h"

#define LOPEZ_SMOOTHING

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

  /* droplets have positive curvature, pay attention to mult */
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

/******************************************************************************/
real VOF::wall_curv_HFnormal_kernel(arr2D & heights,
                                    const arr2D & distances,
                                    const std::vector<real> & dx,
                                    const std::vector<real> & dy,
                                    const bool truedir1, const bool truedir2,
                                    const real mult, const real max_n,
                                    const real cang) {
/***************************************************************************//**
*  \brief Calculate curvature using height-function approach in 3D geometry
*         heights constructed in normal direction.
*
*     output: kappa
*******************************************************************************/
#if 1  
  real & hmm = heights[0][0];
  real & hmc = heights[0][1];
  real & hmp = heights[0][2];
  real & hcm = heights[1][0];
  real & hcc = heights[1][1];
  real & hcp = heights[1][2];
  real & hpm = heights[2][0];
  real & hpc = heights[2][1];
  real & hpp = heights[2][2];

  for(int ii(-hf_set.nof); ii<=hf_set.nof; ++ii)
    for(int jj(-hf_set.nof); jj<=hf_set.nof; ++jj)
      if(heights[ii+hf_set.nof][jj+hf_set.nof]<=0.)
        heights[ii+hf_set.nof][jj+hf_set.nof] = 
          hcc - distances[ii+hf_set.nof][jj+hf_set.nof]*tan(mult*cang);

  /* droplets have positive curvature */
  return calculate_curvature_HF(hmm, hmc, hmp,
                                hcm, hcc, hcp,
                                hpm, hpc, hpp,
                                dx[0], dx[1], dx[2],
                                dy[0], dy[1], dy[2],
                                truedir1, truedir2,
                                mult, max_n); 
#else

#endif
}
