#include "phasechange.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChange::get_M(Scalar & sca) {
/***************************************************************************//**
*  \brief advect mdot in normal direction.
*         call after update (phi (=mdot) is updated)
*         output : sca
*******************************************************************************/
  for_ijk(i,j,k){
    if (phi[i][j][k] == 0.0) {
      sca[i][j][k] = 0.0;
    } else {
      real vol = dV(i,j,k);
       real area = marching_cube(i,j,k);
       if (area!=0.0) {
         sca[i][j][k] = phi[i][j][k]*vol/area;
       } else {
         sca[i][j][k] = 0.0;
       }
    }
  }
  sca.exchange_all();

  return;
}
