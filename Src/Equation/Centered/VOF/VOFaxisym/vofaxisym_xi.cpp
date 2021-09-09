#include "vofaxisym.h"

/***************************************************************************//**
 \brief Calculate the reduced center of mass of given geometry.
        small = mx<mz, negative = nx<0 and vice versa
        source: derived by me (= Lubomir)
*******************************************************************************/
real VOFaxisym::xi_small_pos_triangle(real alpha, real mmx, real mmz) {
  return alpha/3./mmx;
}

real VOFaxisym::xi_small_pos_trapezoid(real alpha, real mmx, real mmz) {
  return (3.*alpha-2.*mmx)/(2.*alpha-mmx)/3.;
}

real VOFaxisym::xi_small_neg_triangle(real alpha, real mmx, real mmz) {
  return (3.*mmx-alpha)/3./mmx;
}

real VOFaxisym::xi_small_neg_trapezoid(real alpha, real mmx, real mmz) {
  return (3.*alpha-mmx)/(2.*alpha-mmx)/3.;
}

real VOFaxisym::xi_large_pos_triangle(real alpha, real mmx, real mmz) {
  return alpha/3./mmx;
}

real VOFaxisym::xi_large_pos_trapezoid(real alpha, real mmx, real mmz) {
  return (2.*alpha-mmz-alpha*(alpha-mmz)/(2.*alpha-mmz))/3./mmx;
}

real VOFaxisym::xi_large_neg_triangle(real alpha, real mmx, real mmz) {
  return (3.*mmx-alpha)/3./mmx;
}

real VOFaxisym::xi_large_neg_trapezoid(real alpha, real mmx, real mmz) {
  return (3.*mmx-2.*alpha+mmz+alpha*(alpha-mmz)/(2.*alpha-mmz))/3./mmx;
}
