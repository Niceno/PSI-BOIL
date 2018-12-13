#include "enthalpytif.h"

/***************************************************************************//**
*  \brief Creates system matrix \f$ [A] \f$.
*
*  It is split into separete parts: discretizaton of innertial and diffusive
*  parts, boundary conditions and. 
*  Finally, it also computes the inverse of the system matrix diagonal.
*******************************************************************************/
void EnthalpyTIF::create_system(const Scalar * diff_eddy) {

   tint_field();


  /*-----------------------------+ 
  |  create system in the core,  |
  |  correct at the boundaries.  |
  +-----------------------------*/
  create_system_innertial();
  create_system_diffusive(diff_eddy);
  create_system_bnd();
#if 0
  std::cout<<"create_system: "<<A.t[1][1][1]<<" "<<A.b[1][1][2]<<"\n";
  exit(0);
#endif

  /*----------------------------------------------+ 
  |  compute invlerse of the central coefficient  |
  +----------------------------------------------*/
  for_ijk(i,j,k) 
    A.ci[i][j][k] = 1.0 / A.c[i][j][k];

  A.ci.exchange();
}

