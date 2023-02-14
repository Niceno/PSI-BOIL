#include "enthalpyfd.h"

/***************************************************************************//**
*  Makes all the necessary preparations for the new time step. That means:
*  fills array \f$ \{f_{old}\} \f$ (fold) with parts of the old innertial term 
*  \f$ \{I\}^{N-1} \f$, old convective term \f$ \{C\}^{N-2} \f$ (cold) and 
*  old diffusive term \f$ \{D\}^{N-1} \f$. 
*
*  Old convective and old diffusive terms are added by calling: 
*  \code convection(&cold) \endcode and \code diffusion() \endcode 
*  respectively.
*
*******************************************************************************/
void EnthalpyFD::new_time_step(const Scalar * diff_eddy) {

  if(conv_ts.Nm1()!=1.0){
    boil::oout<<"enthalpyfd_new_time_step, forward_euler should be used!\n";
    exit(0);
  }
  if(diff_ts.Nm1()!=0.0){
    boil::oout<<"enthalpyfd_new_time_step, backward_euler should be used!\n";
    exit(0);
  }

  convective_time_step();

  /*---------------------------------------+
  |  fold = fold + D                       |
  |  explicit term on the right hand side  |
  +---------------------------------------*/

  if(diff_eddy == NULL){
    if(!laminar){
      boil::oout<<"### enthalpyfd_new_time_step: Error!!!\n";
      boil::oout<<"### mu_t is used for discretization, but is not used for\n";
      boil::oout<<"### new_time_step. Use new_time_step( &mu_t );\n";
      exit(0);
    }
  }
  
  if(diff_ts.Nm1() > 0.0)
    diffusion(diff_eddy);

  return;
}
