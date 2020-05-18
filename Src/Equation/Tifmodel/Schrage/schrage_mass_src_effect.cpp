#include "schrage.h"

/***************************************************************************//**
 *  Incorporates the interfacial resistance effect into the temperature field
 *  mflx is the mass flux value [kg/m^2/s]
 *  mresis is the interfacial resistance to mass transfer, i.e. 1/kin. mobility
******************************************************************************/
void Schrage::mass_src_effect() {
  for_vijk(tif,i,j,k) {
    if(topo->interface(i,j,k)) {
      tif[i][j][k] += mflx[i][j][k]*mresis;
    }
  }

  return;
}
