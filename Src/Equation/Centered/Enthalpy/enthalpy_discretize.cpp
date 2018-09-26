#include "enthalpy.h"

/***************************************************************************//**
*  Discretizes system matrix for the Enthalpy equation.
*******************************************************************************/
void Enthalpy::discretize(const Scalar * eddy) {

  boil::timer.start("enthalpy discretize");

  PropertyDiv * flu_lambda_over_cp = NULL;
  PropertyDiv * sol_lambda_over_cp = NULL;

  flu_lambda_over_cp = new PropertyDiv(flu->lambda(), flu->cp());
  if( solid() )
    sol_lambda_over_cp = new PropertyDiv(sol->lambda(), sol->cp());

  /*-------------------+
  |  no solid defined  |
  +-------------------*/
  if( !solid() ) {
    create_system_innertial(flu->rho());
    create_system_diffusive(flu_lambda_over_cp, NULL, eddy);

  /*-------------------+
  |  solid is defined  |
  +-------------------*/
  } else {
    create_system_innertial(flu->rho(),         sol->rho());
    create_system_diffusive(flu_lambda_over_cp, sol_lambda_over_cp, eddy);
  }

  /*------------------------+ 
  |  correct at boundaries  |
  +------------------------*/  
  create_system_bnd(flu->lambda());

  /*---------------------------+ 
  |  invert the main diagonal  |
  +---------------------------*/  
  for_ijk(i,j,k) 
    A.ci[i][j][k] = 1.0 / A.c[i][j][k];   

  boil::timer.stop("enthalpy discretize");
}

/*-----------------------------------------------------------------------------+
 '$Id: enthalpy_discretize.cpp,v 1.3 2016/02/26 11:07:41 niceno Exp $'/
+-----------------------------------------------------------------------------*/
 
