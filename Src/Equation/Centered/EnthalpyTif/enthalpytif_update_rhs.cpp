#include "enthalpytif.h"

/***************************************************************************//**
*  Called just before solving the linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*  Local array "fnew" represents \f$ \{f\} \f$ in the above equation and
*  contains contribution from time-discretization of various terms (diffusion, 
*  convection, external).
*
*  Note: fext (W) is converted to (W/m3) in this function, taking into account
*  the volume in the framework of finite difference method.
*******************************************************************************/
real EnthalpyTIF::update_rhs() {

  /*--------------------------------------------+ 
  |  update influence from boundary conditions  |
  +--------------------------------------------*/
  phi.bnd_update();

  /*--------------------------+ 
  |  summ all the sources up  |
  +--------------------------*/
  for_ijk(i,j,k) {
    fnew[i][j][k] = fold[i][j][k]
                  + ftif[i][j][k]
                  + cnew[i][j][k] * conv_ts.N()   
                  + fbnd[i][j][k]
                  + fext[i][j][k];
  }

  /*---------------------+ 
  |  this is needed too  |
  +---------------------*/
  if( !solid() ) 
    if(dom->ibody().nccells() > 0) {
      for_ijk(i,j,k)
        if( dom->ibody().off(i,j,k) ) {
          fnew[i][j][k] = phi[i][j][k];
        }
    } /* is there an immersed body */

  return 0.0;
}

