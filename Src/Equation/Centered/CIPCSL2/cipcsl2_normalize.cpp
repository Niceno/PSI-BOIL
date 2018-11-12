#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::normalize(real &nx, real &ny, real &nz) {
/***************************************************************************//**
*  \brief normalize (nx, ny, nz)
*******************************************************************************/

  real nmag = sqrt( nx*nx + ny*ny + nz*nz) + boil::pico;
  nx /= nmag;
  ny /= nmag;
  nz /= nmag;

  return;
}
