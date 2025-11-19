#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Coefficients for the derivatives (to be overwritten by child class).
*******************************************************************************/
real EnthalpyFD::coef_x_m(const real dxm, const real dxp, const real x0) {
  return 2.0/dxm/(dxm+dxp);
}

real EnthalpyFD::coef_x_p(const real dxm, const real dxp, const real x0) {
  return 2.0/dxp/(dxm+dxp);
}

real EnthalpyFD::coef_y_m(const real dxm, const real dxp, const real x0) {
  return 2.0/dxm/(dxm+dxp);
}

real EnthalpyFD::coef_y_p(const real dxm, const real dxp, const real x0) {
  return 2.0/dxp/(dxm+dxp);
}

real EnthalpyFD::coef_z_m(const real dxm, const real dxp, const real x0) {
  return 2.0/dxm/(dxm+dxp);
}

real EnthalpyFD::coef_z_p(const real dxm, const real dxp, const real x0) {
  return 2.0/dxp/(dxm+dxp);
}
