#include "enthalpyfd.h"
#include "def.h"

/***************************************************************************//**
*  Performs a semi-Lagrangian convection step.
*  (convection and diffusion solved separately)
*******************************************************************************/
void EnthalpyFD::convective_time_step() {

  /*---------------------------+
  |  fold = rho * cp * T / dt  |
  +---------------------------*/

  real dti = time->dti();

#ifdef CNEW
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)){
      /* phase change: xor indicates change of phase */
      if(topo->above_interface_old(i,j,k) ^ topo->above_interface(i,j,k)) {
        phi[i][j][k] = tifmodel.Tint(i,j,k);     /* crude code */
      }
    }
  }
  phi.bnd_update();
  phi.exchange();
#endif

  /* no transport in solid */
  if( !solid() ) 
    for_ijk(i,j,k) {
      real c,r;
#ifdef CNEW
      if(topo->above_interface(i,j,k)) {
#else
      if(topo->above_interface_old(i,j,k)) {
#endif
        c = cpl;
      } else {
        c = cpv;
      }
      fold[i][j][k] = c * dV(i,j,k) * phi[i][j][k] * dti;
    }
  /* with transport in solid */
  else {
    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k);
      const real cs = solid()->cp (i,j,k);
      real cf = cpl;
#ifdef CNEW
      if(!topo->above_interface(i,j,k)) {
#else
      if(!topo->above_interface_old(i,j,k)) {
#endif
        cf = cpv;
      }

      fold[i][j][k] = (cf*fV + cs*(1.0-fV)) * dV(i,j,k)
                    * phi[i][j][k] * dti;
    }
  }

  /*-----------------------------------------------------------------------+
  |  fold = fold + C                                                       |
  |  Euler explicit 1st order for convection term                          |
  |  Semi-lagrangian scheme: update convection term, separating diffusion  |
  +-----------------------------------------------------------------------*/
  convection(&cold);
  for_ijk(i,j,k)
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k]; /* conv_ts.Nm1() = 1.5 */

#ifndef CNEW
  /* semi-lagrangian scheme */
  for_ijk(i,j,k){
    if(dom->ibody().on(i,j,k)){
      real c;
      if(topo->above_interface_old(i,j,k)) {
        c = cpl;
      } else {
        c = cpv;
      }
      real t_new = fold[i][j][k] / (c * dV(i,j,k)) / dti;

#if 1
      /* phase change: xor indicates change of phase */
      if(topo->above_interface_old(i,j,k) ^ topo->above_interface(i,j,k)) {
        if( (phi[i][j][k]-tifmodel.Tint(i,j,k))*(t_new-tifmodel.Tint(i,j,k))<=0.0 ){
          t_new = tifmodel.Tint(i,j,k);     /* crude code */
        }
  #if 0
      /* phase does not change */
      } else {
        if( (phi[i][j][k]-tifmodel.Tint(i,j,k))*(t_new-tifmodel.Tint(i,j,k))<0.0 ){
          t_new = tifmodel.Tint(i,j,k);     /* crude code: Is this necessary? */
        }
  #endif
      }
#endif

      phi [i][j][k] = t_new;
      if(topo->above_interface(i,j,k)) {
        c = cpl;
      } else {
        c = cpv;
      }
      fold[i][j][k] = c * dV(i,j,k) * phi[i][j][k] * dti;
    }
  }
  phi.bnd_update();
  phi.exchange();
#endif

  return;
}
