#include "schrage.h"

/***************************************************************************//**
 *  Incorporates the capillary pressure effect into the temperature field
******************************************************************************/
void Schrage::pressure_effect() {
  for_vijk(tif,i,j,k) {
    if(topo->interface(i,j,k)) {
      tif[i][j][k] -= (*dpres)[i][j][k]*tr/rhol/flu->latent(i,j,k);
    }
  }

  return;
}
