#include "centered.h"

/***************************************************************************//**
*  \brief Interface for calling new time step.
*******************************************************************************/
void Centered::new_time_step() {
  if( !solid() ) 
    new_time_step(flu->rho(), NULL);
  else
    new_time_step(flu->rho(), sol->rho());
}

/***************************************************************************//**
*  Makes all the necessary preparations for the new time step. That means:
*  fills array \f$ \{f_{old}\} \f$ (fold) with parts of the old innertial term 
*  \f$ \{I\}^{N-1} \f$, old convective term \f$ \{C\}^{N-2} \f$ (cold) and 
*  old diffusive term \f$ \{D\}^{N-1} \f$. 
*
*  Old convective and old diffusive terms are added by calling: 
*  \code convection(&cold) \endcode and \code diffusion() \endcode 
*  respectivelly.
*******************************************************************************/
void Centered::new_time_step(const Property * f_prop, const Property * s_prop) {

  /*------------------------+
  |      dV  n-1    1  n-2  |
  |  f = -- u    -  - C     |
  |      dt         2       |
  +------------------------*/
  //OPR( conv_ts.Nm1() );
  //OPR( conv_ts.Nm2() );

  /* no transport in solid */
  if( !solid() ) {

    assert(f_prop != NULL);

    for_ijk(i,j,k) {
      const real r = f_prop->value(i,j,k);

      fold[i][j][k] = r * dV(i,j,k) * phi[i][j][k] * time->dti() 
                    + conv_ts.Nm2() * cold[i][j][k]; /* conv_ts.Nm2() = -0.5 */
    }
  }

  /* with transport in solid */
  else {

    assert(f_prop != NULL);
    assert(s_prop != NULL);

    for_ijk(i,j,k) {
      const real rf = f_prop->value(i,j,k);
      const real rs = s_prop->value(i,j,k);
      const real fV = dom->ibody().fV(i,j,k);

      fold[i][j][k] = (rf*fV + rs*(1.0-fV))*dV(i,j,k) 
                    * phi[i][j][k] * time->dti() 
                    + conv_ts.Nm2() * cold[i][j][k]; /* conv_ts.Nm2() = -0.5 */
    }
  }

  /*--------------+
  |       3  n-1  |
  |  f += - C     |
  |       2       |
  +--------------*/
  /* a condition like: if(conv_ts != backward_euler()) would be good */
  convection(&cold, f_prop);
  for_ijk(i,j,k)
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k]; /* conv_ts.Nm1() = 1.5 */

  /*--------------+ 
  |       1  n-1  |
  |  f += - D     |
  |       2       |
  +--------------*/
  /* a condition like: if(diff_ts != backward_euler()) would be good */
  diffusion();

}
