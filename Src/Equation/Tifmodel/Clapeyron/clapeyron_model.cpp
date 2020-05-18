#include "clapeyron.h"

/***************************************************************************//**
 *  Calculates interface temperature field using the Clausius-Clapeyron relation
*******************************************************************************/
void Clapeyron::model() {
  for_vijk(tif,i,j,k) {
    tif[i][j][k] = tr;
    if(topo->interface(i,j,k)) {
      tif[i][j][k] = value(i,j,k);
    }
  }

  return;
}
