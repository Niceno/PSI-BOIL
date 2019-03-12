#include "phasechange.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChange::mdot() {
/***************************************************************************//**
*  \brief advect mdot in normal direction.
*         output : phi
*         tempolary: stmp
*******************************************************************************/

  for_ijk(i,j,k){
    if((iflag[i][j][k] == -1) || (iflag[i][j][k] == 1)){
      if(dom->ibody().on(i,j,k)){
        real mdotc=M[i][j][k];
        real vol = dV(i,j,k);

        /* iso-surface area */
        real area = marching_cube(i,j,k);
        phi[i][j][k] = mdotc * area / vol;
        phi[i][j][k] = mdot_cut(phi[i][j][k],clr[i][j][k]);
      } else {
        phi[i][j][k] = 0.0;
      }
    } else {
      phi[i][j][k] = 0.0;
    }
  }
  phi.exchange_all();

#if 0
  for_i(i) {
    if(fabs(phi[i][1][1])>boil::atto) boil::oout<<"PC::mdot "<<i<<" "<<txv[i][1][1]<<" "<<txl[i][1][1]<<" "<<nx[i][1][1]<<" "<<M[i][1][1]<<" "<<phi[i][1][1]<<boil::endl;
  }
#endif

  return;
}

