#include "concentrationtp.h"
/***************************************************************************//**
*  \brief Creates innertial part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void ConcentrationTP::create_system_innertial() {                     

  /*----------------------------------------------------+ 
  |  compute central coefficient for the system matrix  |
  +----------------------------------------------------*/

  real dti = time->dti();

  for_ijk(i,j,k) {
    real r = rho_dif->value(i,j,k);
    real col_new = (*clr)[i][j][k];

    if(sig==Sign::neg()) col_new = 1.-col_new;
 
    /* gas diffusive innertial */
    A.c[i][j][k] = dV(i,j,k) * dti * r * col_new;
  }
 
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k);
      A.c[i][j][k] *= fV;
    }
  }

  A.c.exchange();
}

