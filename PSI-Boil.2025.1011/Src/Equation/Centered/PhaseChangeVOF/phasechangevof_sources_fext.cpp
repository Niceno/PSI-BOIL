#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::sources_fext() {
/***************************************************************************//**
*  \brief calculate source term for volume conservation
*******************************************************************************/

  for_ijk(i,j,k){
    real mdotc=phi[i][j][k];
    real volc =dV(i,j,k);
    real dt = time->dt();
    fext[i][j][k] = volc/dt*(1.0/rhov-1.0/rhol)*mdotc;
  }
  fext.exchange();

  return;
}
