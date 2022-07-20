#include "enthalpy.h"

/***************************************************************************//**
*  \brief Creates system matrix \f$ [A] \f$.
*
*  It is split into separete parts: discretizaton of innertial and diffusive
*  parts, boundary conditions and. 
*  Finally, it also computes the inverse of the system matrix diagonal.
*******************************************************************************/
void Enthalpy::create_system(const Scalar * diff_eddy) {

  /*-----------------------------+ 
  |  create system in the core,  |
  |  correct at the boundaries.  |
  +-----------------------------*/
  /*-------------------+
  |  no solid defined  |
  +-------------------*/
  if( !solid() ) {
    create_system_innertial(flu->cp());
    create_system_diffusive(flu->lambda(), NULL, diff_eddy);

  /*-------------------+
  |  solid is defined  |
  +-------------------*/
  } else {
    create_system_innertial(flu->cp(),sol->cp());
    create_system_diffusive(flu->lambda(),sol->lambda(), diff_eddy);
  }
  create_system_bnd(flu->lambda());

  /*---------------------------------------------+ 
  |  compute inverse of the central coefficient  |
  +---------------------------------------------*/
  for_ijk(i,j,k) 
    A.ci[i][j][k] = 1.0 / A.c[i][j][k];

  A.ci.exchange();
}

