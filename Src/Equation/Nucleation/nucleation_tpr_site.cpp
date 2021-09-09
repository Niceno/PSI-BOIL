#include "nucleation.h"

/******************************************************************************/
real Nucleation::tpr_site(const int ns ) {
/***************************************************************************//**
*  \brief calculate average temperature of nucleation site
*******************************************************************************/

  real xs = sites[ns].x();
  real ys = sites[ns].y();
  real zs = sites[ns].z();

  real tpr_seed = -boil::unreal;
  if( vf->domain()->contains_xyz( xs, ys, zs) ) {
    /* seed is inside the decomposed domain */
    tpr_seed = (cht->node_tmp_flu())[Comp::k()] /* only k-dir */
                                    [sites[ns].ic()]
                                    [sites[ns].jc()]
                                    [sites[ns].kc()];
  }
  boil::cart.max_real(&tpr_seed);

  return (tpr_seed);
}
