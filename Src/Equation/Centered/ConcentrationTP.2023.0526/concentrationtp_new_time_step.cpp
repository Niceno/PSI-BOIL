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

  if(use_heaviside_instead_of_vf&&!store_vf) {
    for_ijk(i,j,k) {
      vfold[i][j][k] = vfval(i,j,k);
    } 
    store_vf = true;
  }

  /*------------------------------+
  |  fold = vol * rho * eps / dt  |
  +-------------------------------*/
  for_ijk(i,j,k) {
    real col_new = vfval(i,j,k);
    real col_old = vfvalold(i,j,k);
    //real col_new = std::min(1.0,std::max(0.0,clr[i][j][k]));
    //real col_old = std::min(1.0,std::max(0.0,clrold[i][j][k]));
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
  |  Semi-lagrangian scheme: update convection term, separating diffusion  |
  +-----------------------------------------------------------------------*/
  //boil::oout<<"new_time_step:before conv "<<phi[14][33][10]<<"\n";
  convection(&cold);
  //boil::oout<<"conv_ts.Nm1()="<<conv_ts.Nm1()<<"\n"; exit(0);
  for_ijk(i,j,k) { /* conv_ts.Nm1() = 1.0 */
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k];
  }

  //boil::oout<<"new_time_step:after conv "<<phi[14][33][10]<<" "<<matter_sig<<"\n";

#if 1
  real dti = time->dti();

  /* semi-lagrangian scheme */
  for_ijk(i,j,k) {
    real col_new = vfval(i,j,k);
#if 1
    if(matter_sig==Sign::neg()) {
      col_new = 1.-col_new;
    }
#endif
    if(dom->ibody().on(i,j,k)&&heavi->status(i,j,k)!=-matter_sig
       &&col_new>col_crit) {
      real r = rho_dif->value(i,j,k);

      /* gas diffusive innertial */
      real Ac = dV(i,j,k) * dti * r * col_new;

      if(dom->ibody().nccells() > 0) {
        const real fV = dom->ibody().fV(i,j,k);
        Ac *= fV;
      }

      phi[i][j][k] = fold[i][j][k] / Ac;
#if 1
      if(phi[i][j][k]>1.01){
        std::cout<<"new_time_step "<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k]<<" "<<col_new<<" "<<r<<" "<<vfval(i,j,k)<<"\n";
        //exit(0);
      }
#endif
    }
  }
  exit(0);
  //boil::oout<<"new_time_step:after lagran "<<phi[14][33][10]<<"\n";

  extrapolate();
  //boil::oout<<"new_time_step:after extrapolate "<<phi[14][33][10]<<"\n";

  for_ijk(i,j,k) {
    real col_new = vfval(i,j,k);
    if(matter_sig==Sign::neg()) {
      col_new = 1.-col_new;
    }
    if(dom->ibody().on(i,j,k)&&heavi->status(i,j,k)!=-matter_sig
       &&col_new>col_crit) {
      real r = rho_dif->value(i,j,k);

      /* gas diffusive innertial */
      real Ac = dV(i,j,k) * dti * r * col_new;

      if(dom->ibody().nccells() > 0) {
        const real fV = dom->ibody().fV(i,j,k);
        Ac *= fV;
      }

      fold[i][j][k] = Ac * phi[i][j][k];
    }
  }
#endif

  phi.bnd_update();
  phi.exchange();

  /* store vf */
  if(use_heaviside_instead_of_vf) {
    for_ijk(i,j,k) {
      vfold[i][j][k] = vfval(i,j,k);
    }
  }

  boil::timer.stop("concentrationtp new time step");
}
