#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::sources_vfs() {
/***************************************************************************//**
*  \brief calculate source term for volume fraction
*******************************************************************************/

  for_ijk(i,j,k){
    real mdotc=phi[i][j][k];
    vfs[i][j][k]= -(1.0/rhol)*mdotc*sig;
  }
  vfs.exchange();

  return;
}
