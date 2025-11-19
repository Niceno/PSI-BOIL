#include "tif.h"

/***************************************************************************//**
 *  Incorporates the interfacial resistance effect into the temperature field
 *  mflx is the mass flux value [kg/m^2/s]
 *  mresis is the interfacial resistance to mass transfer, i.e. 1/kin. mobility
******************************************************************************/
void TIF::Mass_src_effect() {
  for_vijk(tif,i,j,k) {
    if(Interface(i,j,k)) {
      //if(j==2&&k==2) boil::oout<<"TIF: "<<i<<" "<<tif[i][j][k]<<" "<<tifold[i][j][k]<<" "<<mflx[i][j][k]*mresis<<boil::endl;
      tif[i][j][k] += mflx[i][j][k]*mresis;
    }
  }
}
