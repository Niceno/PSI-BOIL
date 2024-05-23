#include "nucleation.h"

/******************************************************************************/
real Nucleation::tpr_site(const int ns ) {
/***************************************************************************//**
*  \brief calculate average temperature of nucleation site
*******************************************************************************/

  real xs = sites[ns].x();
  real ys = sites[ns].y();
  real zs = sites[ns].z();

  int ii=sites[ns].ic();
  int jj=sites[ns].jc();
  int kk=sites[ns].kc(); // fluid side (Z>0)

  real tpr_seed = -boil::unreal;
  if( vf->domain()->contains_xyz( xs, ys, zs) ) {
    /* seed is inside the decomposed domain */
    //std::cout<<"tpr_site:i,j,k= "<<ii<<" "<<jj<<" "<<kk<<"\n";
    //tpr_seed = (cht->node_tmp_sol())[Comp::k()] // test: output is same as Lubomir
    //tpr_seed = (cht->node_tmp_flu())[Comp::k()] // only k-dir  //Lubomir
    //                                [sites[ns].ic()]
    //                                [sites[ns].jc()]
    //                                [sites[ns].kc()];
    tpr_seed = cht->tmp()[ii][jj][kk-1];
  }
  boil::cart.max_real(&tpr_seed);
  //std::cout<<"tpr_seed= "<<tpr_seed<<"\n";
  //exit(0);

  return (tpr_seed);
}
