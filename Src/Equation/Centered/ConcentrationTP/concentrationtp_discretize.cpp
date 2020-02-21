#include "concentrationtp.h"

/***************************************************************************//**
*  Discretizes system matrix for the Concentration equation.
*******************************************************************************/
void ConcentrationTP::discretize(const Scalar * diff_eddy) {

  boil::timer.start("concentrationtp discretize");

  /* set flag for extrapolation */
  heavi->evaluate_nodes();

  create_system_innertial();
  create_system_diffusive(diff_eddy);

  /*----------------+ 
  |  immersed body?  |
  +-----------------*/  
  if(dom->ibody().nccells() > 0) {
     for_ijk(i,j,k)
       if( dom->ibody().off(i,j,k) ) {
         A.c[i][j][k]  = 1.0;
         A.w[i][j][k]  = 0.0;
         A.e[i][j][k]  = 0.0;
         A.s[i][j][k]  = 0.0;
         A.n[i][j][k]  = 0.0;
         A.b[i][j][k]  = 0.0;
         A.t[i][j][k]  = 0.0;
         A.ci[i][j][k] = 1.0;

         fold[i][j][k] = 0.0;
       }
  } /* is there an immersed body */

  /*-----------+ 
  |  exchange  |
  +-----------*/  
  A.c.exchange();
  A.w.exchange();
  A.e.exchange();
  A.s.exchange();
  A.n.exchange();
  A.b.exchange();
  A.t.exchange();

  /*------------------------+ 
  |  correct at boundaries  |
  +------------------------*/  
  create_system_bnd();

  /*---------------------------+ 
  |  invert the main diagonal  |
  +---------------------------*/  
  for_ijk(i,j,k) {
    if(heavi->status(i,j,k)==1) {
        A.c[i][j][k]  = 1.0;
        A.w[i][j][k]  = 0.0;
        A.e[i][j][k]  = 0.0;
        A.s[i][j][k]  = 0.0;
        A.n[i][j][k]  = 0.0;
        A.b[i][j][k]  = 0.0;
        A.t[i][j][k]  = 0.0;
        A.ci[i][j][k] = 1.0;
        fold[i][j][k] = phi[i][j][k];
    }
    else {
      A.ci[i][j][k] = 1.0 / A.c[i][j][k];
    }
  }
  A.ci.exchange();

  if(diff_eddy != NULL) laminar=false;

  boil::timer.stop("concentrationtp discretize");
}
