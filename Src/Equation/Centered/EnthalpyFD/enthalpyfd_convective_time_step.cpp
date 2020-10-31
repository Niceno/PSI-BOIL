#include "enthalpyfd.h"
#include "def.h"

void EnthalpyFD::convective_time_step() {
  convective_time_step(phi);
}

/***************************************************************************//**
*  Performs a semi-Lagrangian convection step.
*  (convection and diffusion solved separately)
*******************************************************************************/
void EnthalpyFD::convective_time_step(Scalar & sca) {

  /*---------------------------+
  |  fold = rho * cp * T / dt  |
  +---------------------------*/
#ifdef CNEW
  inertial(sca,Old::no);
#else
  inertial(sca,Old::yes);
#endif

  /*-----------------------------------------------------------------------+
  |  fold = fold + C                                                       |
  |  Euler explicit 1st order for convection term                          |
  |  Semi-lagrangian scheme: update convection term, separating diffusion  |
  +-----------------------------------------------------------------------*/
  real dti = time->dti();
  convection(&cold);
  for_ijk(i,j,k)
    fold[i][j][k] += conv_ts.Nm1() * cold[i][j][k]; /* conv_ts.Nm1() = 1.5 */

#ifndef CNEW
  /* semi-lagrangian scheme */
  for_ijk(i,j,k){
    if(dom->ibody().on(i,j,k)){
      real c;
      if(cht.topo->above_interface_old(i,j,k)) {
        c = cht.cpl(i,j,k);
      } else {
        c = cht.cpv(i,j,k);
      }
      real t_new = fold[i][j][k] / (c * dV(i,j,k)) / dti;

#if 1
      /* phase change: xor indicates change of phase */
      if(cht.topo->above_interface_old(i,j,k) ^ cht.topo->above_interface(i,j,k)) {
        if( (phi[i][j][k]-cht.Tint(i,j,k))*(t_new-cht.Tint(i,j,k))<=0.0 ){
          t_new = cht.Tint(i,j,k);     /* crude code */
        }
  #if 0
      /* phase does not change */
      } else {
        if( (phi[i][j][k]-cht.Tint(i,j,k))*(t_new-cht.Tint(i,j,k))<0.0 ){
          t_new = cht.Tint(i,j,k);     /* crude code: Is this necessary? */
        }
  #endif
      }
#endif

      phi [i][j][k] = t_new;
      if(cht.topo->above_interface(i,j,k)) {
        c = cht.cpl(i,j,k);
      } else {
        c = cht.cpv(i,j,k);
      }
      fold[i][j][k] = c * dV(i,j,k) * phi[i][j][k] * dti;
    }
  }
  phi.bnd_update();
  phi.exchange();
#endif

  return;
}
