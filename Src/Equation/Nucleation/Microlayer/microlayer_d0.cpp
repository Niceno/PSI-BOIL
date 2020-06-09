#include "microlayer.h"

/******************************************************************************/
real Microlayer::d0(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief calculate initial thickness of micro layer
*  crude code: assume k-plane
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"microlayer.d0: "<<boil::cart.iam()<<"\n";
#endif

  real rl=distance_from_site(i,j,k);

  real d_return;
  if(!boil::realistic(rl)) {
    d_return = boil::unreal;
  } else if(rl < rmax) {
    d_return = std::max(slope * std::pow(rl, exp_slope), dmicro_min);
  } else {
    d_return = boil::unreal;
  }

  return d_return;
}
