#include "phasechange.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChange::mdot() {
/***************************************************************************//**
*  \brief advect mdot in normal direction.
*         output : phi
*         tempolary: stmp
*******************************************************************************/

  for_vijk(clr,i,j,k){
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
    } else if (iflag[i][j][k] == -3) {
      if(dom->ibody().on(i,j,k)){
        real mdotc=M[i][j][k];
        real vol = dV(i,j,k);
        phi[i][j][k] = mdotc * clr.dSz(i,j,k) / vol;
        phi[i][j][k] = mdot_cut(phi[i][j][k],clr[i][j][k]);
      } else {
        phi[i][j][k] = 0.0;
      }
    } else {
      phi[i][j][k] = 0.0;
    }
  }

  phi.exchange_all();
  //boil::oout<<"mdot[kg/m3/s]= "<<phi[50][50][13]<<"\n";
  //boil::oout<<"Vcell[m3]= "<<dV(50,50,13)<<"\n";
  

  return;
}
