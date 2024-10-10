#include "enthalpyfdhighnu.h"

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
*  buff:   temporarily store updated temperature 
*******************************************************************************/
void EnthalpyFDhighNu::new_time_step(const Scalar * diff_eddy) {

  boil::timer.start("enthalpyfdhighnu new_time_step");

  if(diff_eddy == NULL){
    if(!laminar){
      boil::oout<<"### enthalpyfdhighnu_new_time_step: Error!!!\n";
      boil::oout<<"### mu_t is used for discretization, but is not used for\n";
      boil::oout<<"### new_time_step. Use new_time_step( &mu_t );\n";
      exit(0);
    }
  }

  if(conv_ts.Nm1()!=1.0){
    boil::oout<<"enthalpyfdhighnu_new_time_step, forward_euler must be used!\n";
    exit(0);
  }
  if(diff_ts.Nm1()!=0.0){
    boil::oout<<"enthalpyfdhighnu_new_time_step, backward_euler must be used!\n";
    exit(0);
  }

  /* initial time step or restart */
  if(!store_clrold){
    boil::oout<<"EnthalpyFDhighNu::new_time_step()  initialize clrold"<<"\n";
    for_aijk(i,j,k){
      clrold[i][j][k] = (*clr)[i][j][k];
    }
    store_clrold = true;
  }

  /*---------------------------+
  |  fold = rho * cp * T / dt  |
  +---------------------------*/
  /* no transport in solid */
  if( !solid() ) 
    for_ijk(i,j,k) {
      real c;
      if ((*clr)[i][j][k]>=clrsurf) { // (*clr) is used for system_inertia
      //if (clrold[i][j][k]>=clrsurf) {
        c = cpl;
      } else {
        c = cpv;
      }
      fold[i][j][k] = c * dV(i,j,k) * phi[i][j][k] * time->dti();
    }
  /* with transport in solid */
  else {
    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k);
      const real cs = solid()->cp (i,j,k);
      real cf = cpl;
      if((*clr)[i][j][k]<clrsurf){ // (*clr) is used for system_inertia
      //if(clrold[i][j][k]<clrsurf){
        cf = cpv;
      }

      fold[i][j][k] = (cf*fV + cs*(1.0-fV)) * dV(i,j,k)
                    * phi[i][j][k] * time->dti();
    }
  }

  /*-----------------------------------------------------------------------+
  |  fold = fold + C                                                       |
  |  Euler explicit 1st order for convection term                          |
  |  Semi-laglangian scheme: update convection term, separating diffusion  |
  +-----------------------------------------------------------------------*/
  convection(&cold);
  for_ijk(i,j,k)
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k]; /* conv_ts.Nm1() = 1.0 */

  /* semi-laglangian scheme */
  for_ijk(i,j,k){
    if(dom->ibody().on(i,j,k)){
      // in fluid
      real c;
      if((*clr)[i][j][k]>=clrsurf){  // synchronize with upper definition
        c = cpl;
      } else {
        c = cpv;
      }
      real t_new = fold[i][j][k] / (c * dV(i,j,k)) * time->dt();

      // extrapolate if the phase of the cell center changed
      if( ((*clr)[i][j][k]-clrsurf)*(clrold[i][j][k]-clrsurf) < 0.0){
        real phi_sum = 0.0;
        int num = 0;
        if( ((*clr)[i][j][k]-clrsurf)*(clrold[i-1][j][k]-clrsurf) > 0.0) {
          phi_sum += phi [i-1][j][k];
          num++;
        }
        if( ((*clr)[i][j][k]-clrsurf)*(clrold[i+1][j][k]-clrsurf) > 0.0) {
          phi_sum += phi [i+1][j][k];
          num++;
        }
        if( ((*clr)[i][j][k]-clrsurf)*(clrold[i][j-1][k]-clrsurf) > 0.0) {
          phi_sum += phi [i][j-1][k];
          num++;
        }
        if( ((*clr)[i][j][k]-clrsurf)*(clrold[i][j+1][k]-clrsurf) > 0.0) {
          phi_sum += phi [i][j+1][k];
          num++;
        }
        if( ((*clr)[i][j][k]-clrsurf)*(clrold[i][j][k-1]-clrsurf) > 0.0 &&
            dom->ibody().on(i,j,k-1)){
          phi_sum += phi [i][j][k-1];
          num++;
        }
        if( ((*clr)[i][j][k]-clrsurf)*(clrold[i][j][k+1]-clrsurf) > 0.0) {
          phi_sum += phi [i][j][k+1];
          num++;
        }

        if(num == 0){
          buff[i][j][k] = t_new;
          std::cout<<"Warning:enthalpyfdhighnu_new_time_step:num==0! proc= "
                   <<boil::cart.iam()<<" i,j,k= "<<i<<" "<<j<<" "<<k<<"\n";
        }else{
          buff[i][j][k] = phi_sum/real(num);
        }

      } else {
        buff[i][j][k] = t_new;
      }

      fold[i][j][k] = c * dV(i,j,k) * buff[i][j][k] * time->dti();
    } else {
      // in solid
      buff[i][j][k] = phi[i][j][k];
    }
  }
  phi = buff;
  phi.bnd_update();
  phi.exchange();

  /* store clrold */
  for_aijk(i,j,k){
    clrold[i][j][k] = (*clr)[i][j][k];
  }

  /*---------------------------------------+
  |  fold = fold + D                       |
  |  explicit term on the right hand side  |
  +---------------------------------------*/
  diffusion_fd(diff_eddy);

  boil::timer.stop("enthalpyfdhighnu new_time_step");
}
