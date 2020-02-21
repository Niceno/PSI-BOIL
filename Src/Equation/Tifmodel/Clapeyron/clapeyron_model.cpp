#include "clapeyron.h"

/***************************************************************************//**
 *  Calculates interface temperature field using the Clausius-Clapeyron relation
*******************************************************************************/
void Clapeyron::model() {
  for_vijk(tif,i,j,k) {
    tif[i][j][k] = tr;
    if(Interface(i,j,k)) {
      tif[i][j][k] = value(i,j,k);
    }
  }

  return;
}
