#include "cavitypressure.h"

/***************************************************************************//**
*  Corrects the discretization near the phasic interface using the cavity
*  approximation.
*******************************************************************************/
void CavityPressure::discretize(const Scalar * diff_eddy) {

  /* step 0: base discretization */
  Pressure::discretize_pressure(diff_eddy);

  boil::timer.start("pressure discretize");

  /* step 1: reset gas cells */
  for_ijk(i,j,k) {
    if(!dom->ibody().off_p(i,j,k)&&in_gas(clr[i][j][k])) {
        A.c[i][j][k]  *= boil::pico;
        A.w[i][j][k]  = 0.0;
        A.e[i][j][k]  = 0.0;
        A.s[i][j][k]  = 0.0;
        A.n[i][j][k]  = 0.0;
        A.b[i][j][k]  = 0.0;
        A.t[i][j][k]  = 0.0;
    }
  }

  /* step 2: discretize near-interface liquid cells, using 
             free-surface cavity approximation */
  for_ijk(i,j,k) {
    if(   !dom->ibody().off_p(i,j,k)
       && abs(iflag[i][j][k])<3 && !in_gas(clr[i][j][k])) {
      diff_matrix(i,j,k);
    }
  }

  /* step 3: correct for boundaries & immersed body */

  boil::timer.stop("pressure discretize");

  return;
}
