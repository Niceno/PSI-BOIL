#include "concentrationtp.h"

/***************************************************************************//**
*  Makes all the necessary preparations for the new time step. That means:
*  fills array \f$ \{f_{old}\} \f$ (fold) with parts of the old innertial term 
*  \f$ \{I\}^{N-1} \f$, old convective term \f$ \{C\}^{N-2} \f$ (cold) and 
*  old diffusive term \f$ \{D\}^{N-1} \f$. 
*
*  Old convective and old diffusive terms are added by calling: 
*  \code convection(&cold) \endcode and \code diffusion() \endcode 
*  respectivelly.
*
*  clrold: color function in previous time step
*******************************************************************************/
void ConcentrationTP::new_time_step(const Scalar * diff_eddy) {

  boil::timer.start("concentrationtp new time step");

  if(diff_eddy == NULL) {
    if(!laminar) {
      boil::oout<<"### concentrationtp_new_time_step: Error!!!\n";
      boil::oout<<"### mu_t is used for discretization, but is not used for\n";
      boil::oout<<"### new_time_step. Use new_time_step( &mu_t );\n";
      exit(0);
    }
  }

  if(conv_ts.Nm1()!=1.0) {
    std::cout<<"concentrationtp_new_time_step, forward_euler should be used!\n";
    exit(0);
  }
  if(diff_ts.Nm1()!=0.0) {
    std::cout<<"concentrationtp_new_time_step, backward_euler should be used!\n";
    exit(0);
  }

  /*------------------------------+
  |  fold = vol * rho * eps / dt  |
  +-------------------------------*/
  for_ijk(i,j,k) {
    real col_new = std::min(1.0,std::max(0.0,clr[i][j][k]));
    real col_old = std::min(1.0,std::max(0.0,clrold[i][j][k]));
    //real col_new = clr[i][j][k];
    //real col_old = clrold[i][j][k];
    real r = rho_dif->value(i,j,k);

    if(matter_sig==Sign::neg()) col_old = 1.-col_old;
 
    fold[i][j][k] = r * dV(i,j,k) * time->dti() *
                    phi[i][j][k] * col_old;
  }

  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k);
      fold[i][j][k] *= fV;
    }
  }

  /*-----------------------------------------------------------------------+
  |  fold = fold + C                                                       |
  |  Euler explicit 1st order for convection term                          |
  +-----------------------------------------------------------------------*/
  convection(&cold);
  for_ijk(i,j,k) { /* conv_ts.Nm1() = 1.0 */
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k];
  }

  phi.bnd_update();
  phi.exchange();

  boil::timer.stop("concentrationtp new time step");
}
