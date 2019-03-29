#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::ancillary() {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  boil::timer.start("cipcsl2 ancillary");

  /* wall values extrapolation */
  update_at_walls();

  /*  calculate free surface position */
  cal_fs();

  /* normal vector */
  distfunc(phi,24); 
  gradphic(dist);

  /* calculate area */
  cal_adens();

  boil::timer.stop("cipcsl2 ancillary");

  return;
}

