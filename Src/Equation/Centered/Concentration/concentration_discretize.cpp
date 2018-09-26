#include "concentration.h"

/***************************************************************************//**
*  Discretizes system matrix for the Concentration equation.
*******************************************************************************/
void Concentration::discretize(const Scalar * eddy) {

  boil::timer.start("concentration discretize");

  create_system_innertial(flu->rho());
  create_system_diffusive(flu->gamma(), NULL, eddy);

  /*------------------------+ 
  |  correct at boundaries  |
  +------------------------*/  
  create_system_bnd();

  /*---------------------------+ 
  |  invert the main diagonal  |
  +---------------------------*/  
  for_ijk(i,j,k) 
    A.ci[i][j][k] = 1.0 / A.c[i][j][k];  

  A.ci.exchange();

  boil::timer.stop("concentration discretize");
}

/*-----------------------------------------------------------------------------+
 '$Id: concentration_discretize.cpp,v 1.3 2014/08/06 08:14:31 sato Exp $'/
+-----------------------------------------------------------------------------*/
 
