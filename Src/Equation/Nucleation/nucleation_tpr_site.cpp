#include "nucleation.h"

/******************************************************************************/
real Nucleation::tpr_site(const int ns ) {
/***************************************************************************//**
*  \brief calculate average temperature of nucleation site
*******************************************************************************/

  real tpr_seed = -boil::unreal;

  if (sites[ns].contain_site()) {
    tpr_seed = cht->tmp()[sites[ns].ic_site()]
                         [sites[ns].jc_site()]
                         [sites[ns].kc_site()]; // in solid
  }
  boil::cart.max_real(&tpr_seed);

  if(!boil::realistic(tpr_seed)) {
    boil::oout<<"nucleation:tpr_site: Error!!!\n";
    exit(0);
  }

  return (tpr_seed);
}
