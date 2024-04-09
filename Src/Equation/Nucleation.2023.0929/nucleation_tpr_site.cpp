#include "nucleation.h"
using namespace std;

/******************************************************************************/
real Nucleation::tpr_site(const int ns ) {
/***************************************************************************//**
*  \brief calculate average temperature of nucleation site
*******************************************************************************/

  real xs = sites[ns].x();
  real ys = sites[ns].y();
  //real zs = sites[ns].z();
  real zs = -boil::nano;

  real tpr_seed = -boil::exa;
  if ( tpr->domain()->contains_xyz( xs, ys, zs) ) {
    // seed is inside the decomposed domain
    //tpr_seed = (*tpr)[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()-1]; // crude code: -1
    tpr_seed = (*tpr)[sites[ns].ic()][sites[ns].jc()][sites[ns].kc()]; // crude code: -1
  }
  boil::cart.max_real(&tpr_seed);

  return (tpr_seed);

#if 0
  real atpr=0.0;
  int ntpr=0;

  if ( clr->domain()->contains_xyz( xs, ys, zs) ) {
    if ( clr->domain()->ibody().ncall()==0 ) {
      boil::oout<<"It is not implemented yet.\n";
      boil::oout<<"Condition: Activation of nucleation is controlled by\n";
      boil::oout<<" temperature without solid region.\n";
      exit(0);
    } else {
      for(int cc=0; cc<clr->domain()->ibody().nccells(); cc++){
        int i,j,k;
        clr->domain()->ibody().ijk(cc,&i,&j,&k);
        // (ux,uy,uz) points liquid to solid 
        // crude code!!!
        real ux=clr->domain()->ibody().nwx(i,j,k);
        real uy=clr->domain()->ibody().nwy(i,j,k);
        real uz=clr->domain()->ibody().nwz(i,j,k);
        Dir dir = Dir::undefined();
        if (fabs(uz)>0.707) {
          dir = Dir::kmin();
          if(uz>0.707){
            dir = Dir::kmax();
           }
        } else if (fabs(ux)>0.707) {
          dir = Dir::imin();
          if(ux>0.707){
            dir = Dir::imax();
          }
        } else if (fabs(uy)>0.707) {
          dir = Dir::jmin();
          if(uy>0.707){
            dir = Dir::jmax();
          }
        } else {
          std::cout<<"phasechange_micro: Underdevelopment!!!\n";
          exit(0);
        }

        int iof=0, jof=0, kof=0;
        if(dir == Dir::imin()) iof--; if(dir == Dir::imax()) iof++;
        if(dir == Dir::jmin()) jof--; if(dir == Dir::jmax()) jof++;
        if(dir == Dir::kmin()) kof--; if(dir == Dir::kmax()) kof++;
        // (i ,j ,k ) is in liquid domain
        // (ii,jj,kk) is in solid domain
        int ii=i+iof;
        int jj=j+jof;
        int kk=k+kof;
        real xc = clr->xc(i);
        real yc = clr->yc(j);
        real zc = clr->zc(k);
        real d = sqrt ( pow(xc-xs,2.0) + pow(yc-ys,2.0) + pow(zc-zs,2.0));
        if ( d < rseed ) {
          atpr += (*tpr)[ii][jj][kk];
          ntpr ++;
          //std::cout<<"tpr_site:"<<ii<<" "<<jj<<" "<<kk<<" "
          //         <<(*tpr)[ii][jj][kk]<<"\n";
        }
      }
    }
  }

  boil::cart.sum_real(&atpr);
  boil::cart.sum_int(&ntpr);

  if (ntpr!=0) {
    return(atpr/real(ntpr));
  } else {
    boil::oout<<"tpr_site:Nucleation site doesn't touch to the wall!\n";
    exit(0);
  }
#endif
}
