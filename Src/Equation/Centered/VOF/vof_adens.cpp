#include "vof.h"

/******************************************************************************/
void VOF::cal_adens() {
/***************************************************************************//**
*  \brief Calculate area density in cell.
*         Results: adens
*******************************************************************************/

  //boil::timer.start("vof adens");

  /* cell centered */
  for_vijk(adens,i,j,k) {
    real area = marching_cube_area(i,j,k);
    real volume = dV(i,j,k);

    adens[i][j][k] = area/volume;
  }

#if 0
  real sum(0.0);
  int count(0);
  for_vijk(adens,i,j,k) {
    real sumplus =  adens[i][j][k]*adens.dV(i,j,k);
    if(sumplus>boil::atto) count++;
    sum += sumplus;
  }
  boil::oout<<"VOF::adens "<<count<<" "<<sum<<boil::endl;
#endif

  //boil::timer.stop("vof adens");

  return;
}
