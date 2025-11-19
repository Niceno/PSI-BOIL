#include "enthalpy.h"

/***************************************************************************//**
*  \brief Creates innertial part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void Enthalpy::create_system_innertial(const Property * f_prop,
                                       const Property * s_prop) {

  /*----------------------------------------------------+ 
  |  compute central coefficient for the system matrix  |
  +----------------------------------------------------*/

  real dti = time->dti();

  /* no conduction in solid */
  if( !solid() ) {

    for_ijk(i,j,k) 
      A.c[i][j][k] = dV(i,j,k) * dti * f_prop->value(i,j,k);
 
  /* conduction in solid */
  } else {

    assert(f_prop != NULL);
    assert(s_prop != NULL);

    for_ijk(i,j,k) {
      real fV = dom->ibody().fV(i,j,k); /* fraction in fluid */
      A.c[i][j][k] = dV(i,j,k) * dti * 
                   ( f_prop->value(i,j,k) * fV +
                     s_prop->value(i,j,k) * (1.0-fV) );
    }
  }

  A.c.exchange();
}
