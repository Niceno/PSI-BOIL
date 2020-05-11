#include "antoine.h"

/***************************************************************************//**
 *  Calculates interface temperature field using the Antoine relation
*******************************************************************************/
void Antoine::model() {
  for_vijk(tif,i,j,k) {
    tif[i][j][k] = tr;
    if(Interface(i,j,k)) {
      tif[i][j][k] = value(i,j,k);
    }
  }

  return;
}
