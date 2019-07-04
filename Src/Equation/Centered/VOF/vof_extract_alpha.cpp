#include "vof.h"

/******************************************************************************/
void VOF::extract_alpha() {
/***************************************************************************//**
 \brief Calculate value of alpha in cells
    if there is no interface in the cell, unreal=yotta (=1e+24) is stored.
    plane: vm1*x + vm2*y + vm3*z = alpha
    output: nalpha
    EDIT: I think the alpha is valid in standardized space
*******************************************************************************/

  /* avoid singular cases */
  for_avijk(phi,i,j,k) {
    if(phi[i][j][k]==0.5) phi[i][j][k] += boil::pico;
  }

  /* calculate alpha value in the domain */
  /* assumes positive normal vector in normalized space */
  for_vijk(nalpha,i,j,k) {
    nalpha[i][j][k] = alpha_val(phi[i][j][k],
                                nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  }

  return;
}

/***********************
 * ancillary function
 ***********************/
real VOF::alpha_val(const real c, const real nnx, const real nny, const real nnz) {

  /* degenerate case I */
  if(c<boil::pico||c-1.0>-boil::pico) {
    return boil::unreal;
  }

  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -nnx;
  real vn2 = -nny;
  real vn3 = -nnz;

  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real denom = vm1+vm2+vm3;
  /* degenerate case II */
  if(denom<boil::pico)
    return boil::unreal;

  real qa = 1.0/denom;
  vm1 *= qa;
  vm2 *= qa;
  vm3 *= qa;

  return denom*calc_alpha(c, vm1, vm2, vm3);
}

