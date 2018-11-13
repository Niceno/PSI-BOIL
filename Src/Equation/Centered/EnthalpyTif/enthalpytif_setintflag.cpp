#include "enthalpytif.h"

/***************************************************************************//**
*  \brief Set flag for interface
*  intflag = 1 interface exists in cell
*  intflag = 0 no interface in cell
*******************************************************************************/
void EnthalpyTIF::set_intflag() {

  for_aijk(i,j,k) {
    intflag[i][j][k] = Interface1D_x(i,j,k)
                  || Interface1D_y(i,j,k)
                  || Interface1D_z(i,j,k);
  }

  intflag.exchange();
}	

