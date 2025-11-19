#include "phasechange.h"

/******************************************************************************/
void PhaseChange::sources_clrs() {
/***************************************************************************//**
*  \brief calculate source term for color function
*******************************************************************************/

  for_ijk(i,j,k){
    real mdotc=phi[i][j][k];
    clrs[i][j][k]= -(1.0/rhol)*mdotc;
  }
  clrs.exchange();

  return;
}
