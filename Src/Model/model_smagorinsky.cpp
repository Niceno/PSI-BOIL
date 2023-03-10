#include "model.h"

/***************************************************************************//**
*  Computes eddy viscosity using Smagorinsky SGS model. Unit is [kg/(ms)]. 
*
*  optional (only when damping near wall is necessary): dist_max, dist, y_plus
*  input  c_s: Smagorinsky constant [-]
*         dist_max: maximum distance where damping effect [m]
*         dist: distance function [m]
*  output mu_t: turbulent eddy viscosity [Pa.s] = [kg/(ms)]
*         y_plus: y_plus computed [-]
*******************************************************************************/
void Model::smagorinsky( const Momentum * mom, Scalar * mu_t, real c_s,
             const real *dist_max, const Scalar * dist, Scalar * y_plus) const {

  /* take handy aliases */
  const Vector & uvw   =   mom->val();
  const Matter & fluid = * mom->fluid();
  const Body   & ibody =   mom->domain()->ibody(); 

  assert(mu_t->domain() == uvw.domain() );
  if (y_plus) *y_plus=0.0;
  real damp_sum = 0.0;
  int  damp_n = 0;

  for_vijk( (*mu_t), i, j, k) {
    /*----------------------------------------------+
    |  shear = sqrt( 2 * s_ij * s_ij )       [1/s]  |
    |  s_ij = 1/2 ( du_i/dx_j + du_j/dx_i )  [1/s]  |
    +----------------------------------------------*/
    real shear = uvw.shear_magn(i,j,k);

    /*------------------------------------+
    |  delta = (dx * dy * dz)^(1/3)  [m]  |
    +------------------------------------*/
    real delta = pow((mu_t->dxc(i)*mu_t->dyc(j)*mu_t->dzc(k)),1.0/3.0); 

    /*-----------------------------------------------------------------------+
    |  damping function near wall: Cs = Cs * sqrt(1.0 - exp(-(yplus/25)^3))  |
    |  Eq. (6) in https://doi.org/10.1016/B978-008044114-6/50024-7           |
    +-----------------------------------------------------------------------*/
    real damp = 1.0;
    if (dist) {
      if( ibody.on(i,j,k) ) {
        if ((*dist)[i][j][k] < *dist_max) {
          //std::cout<<"smagorinsky:i,j,k= "<<i<<" "<<j<<" "<<k<<" " 
          //         <<(*dist)[i][j][k]<<" "<<damp_n<<"\n";
          real tx, ty, tz, tau_w, y_pl;
          wall_function( uvw, fluid, * dist, & tau_w, & y_pl,
                         i, j, k, &tx, &ty, &tz );
          (*y_plus)[i][j][k] = y_pl;
          damp = sqrt(1.0 - exp(-pow(y_pl/25.0,3.0)));
          damp_sum += damp;
          damp_n++;
        }
      }
    }

    /*-----------------------------------------------+
    |  mu_t = rho * (Cs * delta)^2 * |s|  [kg/(ms)]  |
    +-----------------------------------------------*/
    real cs_damp = c_s * damp;
    (*mu_t)[i][j][k] = fluid.rho(i,j,k) * pow((cs_damp*delta),2.0)*shear; 
  }

  boil::cart.sum_real(&damp_sum);
  boil::cart.sum_int(&damp_n);
  if (damp_n!=0) {
    boil::oout<<"model_smagorinsky:damp_ave= "<<damp_sum/damp_n<<"\n";
  }

  /*----------------------------------+ 
  |  exchange buffers. important to   |
  |  exchange in the corners as well  |
  +----------------------------------*/
  mu_t->exchange_all();
}
