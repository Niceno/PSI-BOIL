#include "vof.h"

/******************************************************************************/
void VOF::normalize(real &nx, real &ny, real &nz) {
/***************************************************************************//**
*  \brief normalize (nx, ny, nz)
*******************************************************************************/

  real nmag = sqrt( nx*nx + ny*ny + nz*nz) + boil::pico;
  nx /= nmag;
  ny /= nmag;
  nz /= nmag;

  return;
}

/******************************************************************************/
void VOF::normalize_l1(real & nx_l1, real & ny_l1, real & nz_l1,
                       const real nx, const real ny, const real nz) {
/***************************************************************************//**
*  \brief normalize (nx, ny, nz) according to l1 norm, with output to first 3
*******************************************************************************/

  real nmag = sqrt(nx*nx) + sqrt(ny*ny) + sqrt(nz*nz) + boil::pico;
  nx_l1 = nx/nmag;
  ny_l1 = ny/nmag;
  nz_l1 = nz/nmag;

  return;
}
