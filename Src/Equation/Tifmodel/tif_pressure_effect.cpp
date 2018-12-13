#include "tif.h"

/***************************************************************************//**
 *  Incorporates the capillary pressure effect into the temperature field
******************************************************************************/
void TIF::Pressure_effect() {
  for_vijk(tif,i,j,k)
    if(Interface(i,j,k)) {
      tif[i][j][k] -= (*pres)[i][j][k]*tr/rhol/latent;
#if 0
      boil::oout << "ETIF: "<<i<<" "<<(*pres)[i][j][k]<<" "<<(*pres)[i][j][k]*tr/rhol/latent<<boil::endl;
#endif
    }
}
