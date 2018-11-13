#include "phasechangevof.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChangeVOF::mdot() {
/***************************************************************************//**
*         output : phi [kg/m3s]
*******************************************************************************/

#if 0
  boil::oout<<"PCVOF::mdot"<<boil::endl;
#endif

  for_ijk(i,j,k){
   //if(Interface(i,j,k)){
      if(dom->ibody().on(i,j,k)){
        real mdotc=M[i][j][k];

        /* iso-surface area */
        phi[i][j][k] = mdotc * gradclr[i][j][k];
        phi[i][j][k] = mdot_cut(phi[i][j][k],clr[i][j][k]);
      } else {
        phi[i][j][k] = 0.0;
      }
    //} else {
    //  phi[i][j][k] = 0.0;
    //}
  }
  phi.exchange_all();

#if 1
  for_i(i) {
    if(fabs(phi[i][2][2])>boil::atto) boil::oout<<"PCV::mdot "<<i<<" "<<txv[i][2][2]<<" "<<txl[i][2][2]<<" "<<nx[i][2][2]<<" "<<M[i][2][2]<<" "<<phi[i][2][2]<<" "<<gradclr[i][2][2]<<" "<<tpr[i][2][2]<<" "<<tpr[i+1][2][2]<<" "<<tpr[i+2][2][2]<<boil::endl;
  }
#endif

  return;
}

