#include "model.h"

/***************************************************************************//**
*  Werner and Wengle Model
*  u+ = y+         if y+ < 11.84
*  u+ = A(y+)^B    otherwise
*
*  Note: u+ = u / u_tau       --┐
*        y+ = y * u_tau / nu  ----> calculate u_tau
*******************************************************************************/
void Model::wall_function( const Vector & uvw, const Matter & fluid,
                           const Scalar & dist, 
                           real * tau_w, real * y_pl,
                           const int i, const int j, const int k,
                           real * tx, real * ty, real * tz ) const {

  /*---------------------------+
  |  constants from power law  |
  +---------------------------*/
  const real a_pow = 8.3;
  const real b_pow = 1.0/7.0;

  /* get normal to the wall */
  real nx, ny, nz; dist.grad_abs(i, j, k, &nx, &ny, &nz);

  /* estimate tangential velocity component */
  real ux, uy, uz; uvw.central( i, j, k, &ux, &uy, &uz );
  const real ut = tangential_velocity(ux, uy, uz, nx, ny, nz, tx, ty, tz);

  /* get distance to the wall */
  const real d  = dist[i][j][k];
 
  /* and kinematic viscosity */
  const real nu = fluid.mu(i,j,k) / fluid.rho(i,j,k);

  /* try linear */
  real y_plus = sqrt( ut * d / nu );

  /* low re */
  if( y_plus < 11.84 ) {
    *y_pl  = y_plus;
    *tau_w = fluid.rho(i,j,k) * ut*nu/d;

  /* high re */
  } else {
    y_plus = pow( ut*d/a_pow/nu, 1.0/(1.0+b_pow) );
    *y_pl  = y_plus;
    real tmp = ut/a_pow * pow( nu/d , b_pow );
    real u_tau = pow (tmp, 1.0/(1.0+b_pow));
    *tau_w = fluid.rho(i,j,k) * pow( u_tau, 2.0 ); 
    //*tau_w = fluid.rho(i,j,k) * pow( tmp, (2.0/(1.0+b_pow) ) ); 
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: model_wall_function.cpp,v 1.7 2015/08/07 07:33:33 sato Exp $'/
+-----------------------------------------------------------------------------*/
