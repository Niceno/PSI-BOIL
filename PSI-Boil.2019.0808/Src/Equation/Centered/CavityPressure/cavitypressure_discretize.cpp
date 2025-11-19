#include "cavitypressure.h"

/***************************************************************************//**
*  Corrects the discretization near the phasic interface using the cavity
*  approximation.
*******************************************************************************/
void CavityPressure::discretize(const Scalar * diff_eddy) {

  boil::timer.start("cavitypressure discretize");

  /* initialize explicit term */
  fold = 0.0;

  /* step 1: discretize cells, using 
             free-surface cavity approximation */
  for_ijk(i,j,k) {
    diff_matrix(i,j,k);
  }

  /* step 2: reset gas and immersed body cells */
  for_ijk(i,j,k) {
    if(dom->ibody().off_p(i,j,k)||in_gas(i,j,k)) {
        A.c[i][j][k]  *= boil::pico;
        A.w[i][j][k]  = 0.0;
        A.e[i][j][k]  = 0.0;
        A.s[i][j][k]  = 0.0;
        A.n[i][j][k]  = 0.0;
        A.b[i][j][k]  = 0.0;
        A.t[i][j][k]  = 0.0;
        fold[i][j][k] = 0.0;
    }
  }

  /* step 3: correct for boundaries & immersed body */
  create_system_bnd();

  /* step 4: calculate inverse central coefficient */
  for_ijk(i,j,k) A.ci[i][j][k] = 1.0 / A.c[i][j][k];

  boil::timer.stop("cavitypressure discretize");

  return;
}
