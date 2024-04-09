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
   if(Interface(i,j,k)){
      if(dom->ibody().on(i,j,k)){
        real mdotc=M[i][j][k];

        /* iso-surface area */
        phi[i][j][k] = mdotc * adens[i][j][k];
        phi[i][j][k] = mdot_cut(phi[i][j][k],clr[i][j][k]);
      } else {
        phi[i][j][k] = 0.0;
      }
    } else {
      phi[i][j][k] = 0.0;
    }
  }
#if 1
  phi.bnd_update();
  phi.exchange_all();
#else
  /* cut */
  mdot_cut();
#endif

  return;
}

