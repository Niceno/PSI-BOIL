#include "enthalpytif.h"

/***************************************************************************//**
 *  Incorporates the capillary pressure effect into the temperature field
******************************************************************************/
void EnthalpyTIF::Pressure_effect(const Scalar & heaviside) {
  for_ijk(i,j,k)
    if(Interface(i,j,k,heaviside)) {
      tif[i][j][k] -= (*pres)[i][j][k]*tsat/rhol/latent;
#if 0
      boil::oout << "ETIF: "<<i<<" "<<(*pres)[i][j][k]<<" "<<(*pres)[i][j][k]*tsat/rhol/latent<<boil::endl;
#endif
    }
}
