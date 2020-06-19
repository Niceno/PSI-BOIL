#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::sources_vfs() {
/***************************************************************************//**
*  \brief calculate source term for volume fraction
*******************************************************************************/

  for_ijk(i,j,k){
    real vfsc = phi[i][j][k];
    vfsc *= -(1.0/rhol)*matter_sig;
    vfs[i][j][k]= vfs_cut(vfsc,vf[i][j][k]);
  }
  vfs.exchange();

  return;
}
