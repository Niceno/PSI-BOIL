#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Creates innertial part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void EnthalpyFD::create_system_innertial() {                     

  /*----------------------------------------------------+ 
  |  compute central coefficient for the system matrix  |
  +----------------------------------------------------*/

  real dti = time->dti();

  /* no conduction in solid */
  if( !solid() ) {

    for_ijk(i,j,k){
      real r,c;
      if((*clr)[i][j][k]>=clrsurf){
        c = cpl;
      } else {
        c = cpv;
      }
      A.c[i][j][k] = dV(i,j,k) * dti * c;
    }
 
  /* conduction in solid */
  } else {

    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k); /* fraction in fluid */
      real c=cpl;
      if((*clr)[i][j][k]<clrsurf){
        c = cpv;
      }

      A.c[i][j][k] = dV(i,j,k) * dti * 
                   ( c * fV + solid()->cp(i,j,k) * (1.0-fV) );
    }
  }

  A.c.exchange();
}
