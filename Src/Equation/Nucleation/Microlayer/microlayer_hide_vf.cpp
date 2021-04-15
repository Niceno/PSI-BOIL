#include "microlayer.h"

/******************************************************************************/
void Microlayer::hide_vf() {
/***************************************************************************//**
*  \brief hide the volume fraction field in ml before advancing the itm
*******************************************************************************/

  for(int cc=0; cc<cht->topo->domain()->ibody().nccells(); cc++){
    int i,j,k;
    cht->topo->domain()->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=cht->topo->domain()->ibody().nwx(i,j,k);
    real uy=cht->topo->domain()->ibody().nwy(i,j,k);
    real uz=cht->topo->domain()->ibody().nwz(i,j,k);
    Comp mcomp;
    Sign sig = Sign::neg();
    Dir d = Dir::undefined();
    if (fabs(uz)>0.707) {
      mcomp = Comp::k();
      d = Dir::kmin();
      if(uz>0.707){
        d = Dir::kmax();
        sig = Sign::pos();
      }
    } else if (fabs(ux)>0.707) {
      mcomp = Comp::i();
      d = Dir::imin();
      if(ux>0.707){
        d = Dir::imax();
        sig = Sign::pos();
      }
    } else if (fabs(uy)>0.707) {
      mcomp = Comp::j();
      d = Dir::jmin();
      if(uy>0.707){
        d = Dir::jmax();
        sig = Sign::pos();
      }
    } else {
      boil::oout<<"microlayer_area_effect: Underdevelopment!!!\n";
      exit(0);
    }

    /* only cells with color < 0.5 */
    if(in_vapor(i,j,k)) {

      /* neighbouring cells */
      std::array< std::array<int,3>, 4> nbrs
        = {{ {{i,j,k}}, {{i,j,k}}, {{i,j,k}}, {{i,j,k}} }};
      /*  ^--- double braces necessary due to the way C++ standard works */

      /* indices in directions perpendicular to mcomp */
      int jdx(1), kdx(2);
      if(mcomp==Comp::i()) {
      } else if(mcomp==Comp::j()) {
        jdx = 2;
        kdx = 0;
      } else {
        jdx = 0;
        kdx = 1;
      }

      /* 0,-1,0 */
      nbrs[0][jdx] += -1;
      /* 0,+1,0 */
      nbrs[1][jdx] +=  1;
      /* 0,0,-1 */
      nbrs[2][kdx] += -1;
      /* 0,0,+1 */
      nbrs[3][kdx] +=  1;

      /* exclude cells, which have liquid as neighbours */
      bool cond = true;
      for(auto & nc : nbrs) {
        //boil::oout<<"HERE "<<i<<" "<<j<<" "<<k<<" "<<nc[0]<<" "<<nc[1]<<" "<<nc[2]<<boil::endl;
        cond &= in_vapor(nc[0],nc[1],nc[2]);
      }
      //boil::oout<<cond<<boil::endl;
      if(!cond)
        continue;

      /* hide microlayer */
      if(boil::realistic(dmicro[i][j][k])&&dmicro[i][j][k]>0.) {
        (*vf)[i][j][k] = 0.0;
        (*cht->topo->clr)[i][j][k] = 0.0;
        cht->topo->vfold[i][j][k] = 0.0;
      }

    }/* in vapor */
  } /* ibody cells */

  return;
}
