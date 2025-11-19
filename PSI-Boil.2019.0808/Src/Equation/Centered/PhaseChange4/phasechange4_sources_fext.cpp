#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::sources_fext() {
/***************************************************************************//**
*  \brief calculate source term for volume conservation
*******************************************************************************/

  for_ijk(i,j,k){
    real mdotc=phi[i][j][k];
    real volc =dV(i,j,k);
    real dt = time->dt();
    fext[i][j][k] = volc/dt*(1.0/cht.rhov(i,j,k)-1.0/cht.rhol(i,j,k))
                   *mdotc*matter_sig;
  }
  fext.exchange();

  return;
}
