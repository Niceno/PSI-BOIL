#include "enthalpytif.h"

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
void EnthalpyTIF::new_time_step(const Scalar * diff_eddy) {

  if(diff_eddy == NULL){
    if(!laminar){
      boil::oout<<"### enthalpytif_new_time_step: Error!!!\n";
      boil::oout<<"### mu_t is used for discretization, but is not used for\n";
      boil::oout<<"### new_time_step. Use new_time_step( &mu_t );\n";
      exit(0);
    }
  }

  if(conv_ts.Nm1()!=1.0){
    std::cout<<"enthalpytif_new_time_step, forward_euler should be used!\n";
    exit(0);
  }
  if(diff_ts.Nm1()!=0.0){
    std::cout<<"enthalpytif_new_time_step, backward_euler should be used!\n";
    exit(0);
  }

  /* initial time step or restart */
  if(!store_clrold){
    boil::oout<<"EnthalpyFD::new_time_step()  initialize clrold"<<"\n";
    for_aijk(i,j,k){
      clrold[i][j][k] = (*clr)[i][j][k];
    }
    if(fs) {
      for_m(m)
        for_avmijk(fsold,m,i,j,k)
          fsold[m][i][j][k] = (*fs)[m][i][j][k];
    }
    store_clrold = true;
  }

  /*---------------------------+
  |  fold = rho * cp * T / dt  |
  +---------------------------*/
  /* no transport in solid */
  if( !solid() ) 
    for_ijk(i,j,k) {
      real c,r;
      if (clrold[i][j][k]>=clrsurf) {
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
      if(clrold[i][j][k]<clrsurf){
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
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k]; /* conv_ts.Nm1() = 1.5 */

  /* semi-lagrangian scheme */
  for_ijk(i,j,k){
    if(dom->ibody().on(i,j,k)){
      real c;
      if(clrold[i][j][k]>=clrsurf){
        c = cpl;
      } else {
        c = cpv;
      }
      real t_new = fold[i][j][k] / (c * dV(i,j,k)) * time->dt();

#if 0
      // phase change
      if( ((*clr)[i][j][k]-clrsurf)*(clrold[i][j][k]-clrsurf) < 0.0){
        if( (phi[i][j][k]-tsat)*(t_new-tsat)<=0.0 ){
          t_new = tsat;     /* crude code */
        }
      // phase does not change
      } else {
        if( (phi[i][j][k]-tsat)*(t_new-tsat)<0.0 ){
          t_new = tsat;     /* crude code: Is this necessary? */
        }
      }
#endif

      phi [i][j][k] = t_new;
      if((*clr)[i][j][k]>=clrsurf){
        c = cpl;
      } else {
        c = cpv;
      }
      fold[i][j][k] = c * dV(i,j,k) * phi[i][j][k] * time->dti();
    }
  }
  phi.bnd_update();
  phi.exchange();

  /* store clrold */
  for_aijk(i,j,k){
    clrold[i][j][k] = (*clr)[i][j][k];
  }
  if(fs) {
    for_m(m)
      for_avmijk(fsold,m,i,j,k)
        fsold[m][i][j][k] = (*fs)[m][i][j][k];
  }

  /*---------------------------------------+
  |  fold = fold + D                       |
  |  explicit term on the right hand side  |
  +---------------------------------------*/
  
  diffusion_fd(diff_eddy);

}

