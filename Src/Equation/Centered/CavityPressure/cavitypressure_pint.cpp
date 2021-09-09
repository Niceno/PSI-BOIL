#include "cavitypressure.h"

/***************************************************************************//**
 *  Returns interface pressure modelled by the cavity 
*******************************************************************************/
real CavityPressure::Pint(const int i, const int j, const int k) {
  /* in psi-boil, bubbles have negative curvature:
   *           pl = pg + sigma*kappa;
   * multiplication by sig guarantees correct value under inversion. 
   * (this is normally handled by the normal vector inversion) */
  return matter_sig*sigma->value(i,j,k)*(*kappa)[i][j][k]+cavity_pressure;
}

/***************************************************************************//**
 *  Returns pressure modelled by the cavity 
*******************************************************************************/
real CavityPressure::Pcav(const int i, const int j, const int k) {
  return cavity_pressure;
}
