#include "microlayer.h"

/******************************************************************************/
real Microlayer::d0(const int i, const int j, const int k) const {
/***************************************************************************//**
*  \brief calculate initial thickness of micro layer
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"microlayer.d0: "<<boil::cart.iam()<<"\n";
#endif

  real rl=distance_from_site(i,j,k);

  real d_return;
  if(!boil::realistic(rl)) {
    d_return = boil::unreal;
  } else if(rl < rmax) {
    d_return = std::min(dmicro_max,
                        std::max(slope * std::pow(rl, exp_slope), dmicro_min));
  } else {
    d_return = boil::unreal;
  }

#if 0
  if (i==26&&j==21&&k==15)
  std::cout<<"microlayer_d0:26,21,15 "<<rl<<" "<<d_return<<"\n";
#endif
  return d_return;
}

/******************************************************************************/
real Microlayer::d0max(const Comp mcomp, const int i, 
                       const int j, const int k) const {
/***************************************************************************//**
*  \brief compare to vf limit
*******************************************************************************/
  /* the sign for face distance does not matter */
  if(matter_sig==Sign::pos()) {
    return threshold_c*2.*cht->distance_face(Sign::neg(),mcomp,i,j,k);
  } else {
    return (1.-threshold_c)*2.*cht->distance_face(Sign::neg(),mcomp,i,j,k);
  }
}
