#include "electpoten.h"

/***************************************************************************//**
*  Since discretization of pressure-Poisson equation is different from 
*  discretizations of other variables in a fact that diffusion-like
*  neighbouring coefficients depend on the inverse of physical property
*  (density), a specialized function is written here. 
*
*  The only thing useful from the parent is handling of boundary conditions.
*******************************************************************************/
void ElectPoten::discretize(const Scalar * diff_eddy) {

  boil::timer.start("pressure discretize");

  Comp m;
  real a_w, a_e, a_s, a_n, a_b, a_t;
  real sigm, sigp;

  /*----------------------------+
  |  neighbouring coefficients  |
  +----------------------------*/

  /* coefficients in i direction (w and e) */
  m = Comp::u();
  for_ijk(i,j,k) {
    /* linear */
    sigm = fluid()->sigma_e(m,i,  j,k);
    sigp = fluid()->sigma_e(m,i+1,j,k);
    a_w = dSx(Sign::neg(),i,j,k);
    a_e = dSx(Sign::pos(),i,j,k);
    A.w[i][j][k] = a_w / dxw(i) * sigm;
    A.e[i][j][k] = a_e / dxe(i) * sigp;
  }

  /* coefficients in j direction (s and n) */
  m = Comp::v();
  for_ijk(i,j,k) {
    /* linear */
    sigm = fluid()->sigma_e(m,i,j,  k);
    sigp = fluid()->sigma_e(m,i,j+1,k);
    a_s = dSy(Sign::neg(),i,j,k);
    a_n = dSy(Sign::pos(),i,j,k);
    A.s[i][j][k] = a_s / dys(j) * sigm;
    A.n[i][j][k] = a_n / dyn(j) * sigp;
  }

  /* coefficients in k direction (b and t) */
  m = Comp::w();
  for_ijk(i,j,k) {
    /* linear */
    sigm = fluid()->sigma_e(m,i,j,k);
    sigp = fluid()->sigma_e(m,i,j,k+1);
    a_b = dSz(Sign::neg(),i,j,k);
    a_t = dSz(Sign::pos(),i,j,k);
    A.b[i][j][k] = a_b / dzb(k) * sigm;
    A.t[i][j][k] = a_t / dzt(k) * sigp;
  }

  /*----------------------+
  |  central coefficient  |
  +----------------------*/
  /* since solid and fluid cells are blended during coarsening,
     it is necessary for the central coefficient to maintain
     magnitude negligible wrt fluid ones (comment 1/2) */
  for_ijk(i,j,k) A.c[i][j][k] = A.w[i][j][k] + A.e[i][j][k]
                              + A.s[i][j][k] + A.n[i][j][k]
                              + A.b[i][j][k] + A.t[i][j][k];

  /*------------------------+
  |  boundaries correction  |
  +------------------------*/
  /* 19-11-28: moved here to allow other bcs than zero-Neumann */
  Centered::create_system_bnd();

  for_ijk(i,j,k) A.ci[i][j][k] = 1.0 / A.c[i][j][k];

  boil::timer.stop("pressure discretize");
}
