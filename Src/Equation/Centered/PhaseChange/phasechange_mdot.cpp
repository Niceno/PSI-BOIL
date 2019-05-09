#include "phasechange.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChange::mdot() {
/***************************************************************************//**
*  \brief advect mdot in normal direction.
*         output : phi
*         tempolary: stmp
*******************************************************************************/

  heavi.calculate_adens();
  //boil::plot->plot(adens,"adens",time->current_step());
  //exit(0);

  for_ijk(i,j,k){
    if((iflag[i][j][k] == -1) || (iflag[i][j][k] == 1)){
      if(dom->ibody().on(i,j,k)){
        real mdotc=M[i][j][k];
        real vol = dV(i,j,k);

        /* iso-surface area */
        real ardens = adens[i][j][k];
        phi[i][j][k] = mdotc * ardens;
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
  for_ijk(i,j,k) {
    if(fabs(phi[i][j][k])>boil::atto) boil::oout<<"PC::mdot "<<i<<" "<<txv[i][j][k]<<" "<<txl[i][j][k]<<" "<<nx[i][j][k]<<" "<<M[i][j][k]<<" "<<phi[i][j][k]<<" "<<adens[i][j][k]<<boil::endl;
  }
#endif

  return;
}

