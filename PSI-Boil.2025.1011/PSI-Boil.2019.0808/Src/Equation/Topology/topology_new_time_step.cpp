#include "topology.h"

/******************************************************************************/
void Topology::new_time_step() {
/***************************************************************************//**
*  \brief store old topology.
*******************************************************************************/

  for_avijk(clrold,i,j,k) {
      clrold[i][j][k]   = (*clr)[i][j][k];
      vfold[i][j][k]   = (*vf)[i][j][k];
      iflagold[i][j][k] = (*iflag)[i][j][k];
  }
  for_m(m)
    for_avmijk(fsold,m,i,j,k)
      fsold[m][i][j][k] = (*fs)[m][i][j][k];
  
  return;
}
