#include "phasechange.h"

/******************************************************************************/
void PhaseChange::sources_sum() {
/***************************************************************************//**
*  \brief calculate source terms
*******************************************************************************/

  // sum positive-mdot and negative-mdot
  // mdot: positive for boiling, negative for condensation
  smdot_pos = 0.0;
  smdot_neg = 0.0;
  for_vijk(phi,i,j,k){
#ifdef IB
    if(dom->ibody().off(i,j,k))continue;
#endif
    real volc =dV(i,j,k);
    if(phi[i][j][k]>0.0) smdot_pos += phi[i][j][k]*volc;
    if(phi[i][j][k]<0.0) smdot_neg += phi[i][j][k]*volc;
  }
  boil::cart.sum_real(&smdot_pos);
  boil::cart.sum_real(&smdot_neg);

  return;
}
