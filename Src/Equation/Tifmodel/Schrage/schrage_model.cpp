#include "schrage.h"

/***************************************************************************//**
 *  Calculates interface temperature field using Schrage model
*******************************************************************************/
void Schrage::model() {
  for_vijk(tif,i,j,k) {
    tif[i][j][k] = tr;
  }

  if(dpres)
    pressure_effect();

  mass_src_effect();

  return;
}
