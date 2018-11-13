#include "enthalpytif.h"

/***************************************************************************//**
 *  Incorporates the interfacial resistance effect into the temperature field
 *  mflx is the mass flux value [kg/m^2/s]
 *  mresis is the interfacial resistance to mass transfer, i.e. 1/kin. mobility
******************************************************************************/
void EnthalpyTIF::Mass_src_effect(const Scalar & heaviside) {
  for_ijk(i,j,k)
    if(Interface(i,j,k,heaviside))
      tif[i][j][k] += (*mflx)[i][j][k]*mresis;
}
