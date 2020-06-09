#include "nucleation.h"

/******************************************************************************/
real Nucleation::clr_site(const int ns ) {
/***************************************************************************//**
*  \brief obtain color at nucleation site
*******************************************************************************/

  real clr_seed = -boil::unreal;
  real xs = sites[ns].x();
  real ys = sites[ns].y();
  real zs = sites[ns].z();
  //std::cout<<"clr_site: "<<sites[ns].ic()<<" "<<sites[ns].jc()<<" "
  //         <<sites[ns].kc()<<" "<<boil::cart.iam()<<"\n";
  if ( clr->domain()->contains_xyz( xs, ys, zs) ) {
    // seed is inside the decomposed domain
    clr_seed = (*clr)[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()];
  }
  boil::cart.max_real(&clr_seed);

  if(!boil::realistic(clr_seed)) {
    boil::oout<<"nucleation:clr_site: Error!!! Cannot find seed point!\n";
    exit(0);
  }

  return clr_seed;
}
