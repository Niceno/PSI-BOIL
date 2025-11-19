#include "vof.h"

/******************************************************************************/
void VOF::cal_adens() {
/***************************************************************************//**
*  \brief Calculate area density in cell.
*         Results: adens
*******************************************************************************/

  //boil::timer.start("vof adens");
  heavi.calculate_adens();
  //boil::timer.stop("vof adens");

  return;
}
