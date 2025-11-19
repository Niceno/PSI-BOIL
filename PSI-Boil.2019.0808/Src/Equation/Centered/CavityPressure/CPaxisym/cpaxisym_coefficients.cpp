#include "cpaxisym.h"

/***************************************************************************//**
*  \brief Coefficients for the derivatives (to be overwritten by child class).
*******************************************************************************/
real CPaxisym::coef_x_m(const real dxm, const real dxp, const real x0) {
  return 1.0/dxm/(dxm+dxp) * (1.0+(x0-dxm)/x0);
}

real CPaxisym::coef_x_p(const real dxm, const real dxp, const real x0) {
  return 1.0/dxp/(dxm+dxp) * (1.0+(x0+dxp)/x0);
}

real CPaxisym::coef_y_m(const real dxm, const real dxp, const real x0) {
  return 0.0;
}

real CPaxisym::coef_y_p(const real dxm, const real dxp, const real x0) {
  return 0.0;
}

real CPaxisym::coef_z_m(const real dxm, const real dxp, const real x0) {
  return 2.0/dxm/(dxm+dxp);
}

real CPaxisym::coef_z_p(const real dxm, const real dxp, const real x0) {
  return 2.0/dxp/(dxm+dxp);
}
