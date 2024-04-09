#include "nucleation.h"
using namespace std;

/******************************************************************************/
real Nucleation::clr_site(const int ns ) {
/***************************************************************************//**
*  \brief calculate average temperature of nucleation site
*******************************************************************************/

  real clr_seed = -boil::exa;
  real xs = sites[ns].x();
  real ys = sites[ns].y();
  //real zs = sites[ns].z();
  real zs = -boil::nano;
  //std::cout<<"clr_site: "<<sites[ns].ic()<<" "<<sites[ns].jc()<<" "
  //         <<sites[ns].kc()<<" "<<boil::cart.iam()<<"\n";
  if ( clr->domain()->contains_xyz( xs, ys, zs) ) {
    // seed is inside the decomposed domain
    clr_seed = (*clr)[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()+1];
  }
  boil::cart.max_real(&clr_seed);

  if(clr_seed<=-boil::exa) {
    boil::oout<<"nucleation:clr_site: Error!!! Cannot find seed point!\n";
    exit(0);
  }

  return clr_seed;
}
