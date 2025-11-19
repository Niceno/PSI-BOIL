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
      if(cht.above_interface(i,j,k,Old::no)) {
        c = cht.cpl(i,j,k);
      } else {
        c = cht.cpv(i,j,k);
      }
      A.c[i][j][k] = dV(i,j,k) * dti * c;
    }
 
  /* conduction in solid */
  } else {

    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k); /* fraction in fluid */
      real c=cht.cpl(i,j,k);
      if(!cht.above_interface(i,j,k,Old::no)) {
        c = cht.cpv(i,j,k);
      }

      A.c[i][j][k] = dV(i,j,k) * dti * 
                   ( c * fV + solid()->cp(i,j,k) * (1.0-fV) );
    }
  }

  A.c.exchange();
}

