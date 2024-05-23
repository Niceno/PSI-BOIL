#include "nucleation.h"

/******************************************************************************/
real Nucleation::clr_site(const int ns ) {
/***************************************************************************//**
*  \brief obtain color at nucleation site
*******************************************************************************/

  real clr_seed = -boil::unreal;
  if (sites[ns].contain_site()) {
    clr_seed = (*(cht->topo->clr))[sites[ns].ic_site()]
                                  [sites[ns].jc_site()]
                                  [sites[ns].kc_site()+1];  // fluid phase
  }
  boil::cart.max_real(&clr_seed);

  if(!boil::realistic(clr_seed)) {
    boil::oout<<"nucleation:clr_site: Error!!! Cannot find seed point!\n";
    exit(0);
  }

  return clr_seed;
}
