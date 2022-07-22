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
#if 1
    tpr_seed = (cht->node_tmp_flu())[Comp::k()] /* only k-dir */
                                    [sites[ns].ic()]
                                    [sites[ns].jc()]
                                    [sites[ns].kc()];
#else
    std::cout<<"tpr_site1: i,j,k=" <<(cht->tmp())[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()]
       <<" "<<sites[ns].ic()<<" "<<sites[ns].jc()<<" "<<sites[ns].kc()<<" "
       <<sites[ns].z()<<"\n";
    std::cout<<"tpr_site2: "<<(cht->tmp())[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()-1]
	     <<" "<<(cht->node_tmp_flu())[Comp::k()][sites[ns].ic()][sites[ns].jc()][sites[ns].kc()]
	     <<" "<<(cht->node_tmp_flu())[Comp::k()][sites[ns].ic()][sites[ns].jc()][sites[ns].kc()-1]
	     <<" "<<(cht->node_tmp_flu())[Comp::k()][sites[ns].ic()][sites[ns].jc()][sites[ns].kc()+1]
	     <<" "<<(cht->node_tmp_flu()).zc(Comp::k(),sites[ns].kc())<<"\n";
    tpr_seed = (cht->tmp())[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()-1]; //crude code
#endif
  }
  boil::cart.max_real(&tpr_seed);
  //std::cout<<"tpr_seed= "<<tpr_seed<<"\n";

  return (tpr_seed);
}
