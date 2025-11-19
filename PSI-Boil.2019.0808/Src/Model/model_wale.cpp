#include "model.h"

/***************************************************************************//**
*  Computes eddy viscosity using WALE SGS model. Unit is [kg/(ms)].
*
*  F.Nicoud and F.Ducros, Subgrid-scale stress..., Flow Turbulence and 
*  Combustion, 62, 183-200, 1999
*
*******************************************************************************/
void Model::wale( const Momentum * mom, Scalar * mu_t, real c_w ) const {

  /* take handy aliases */
  const Vector & uvw   =   mom->val();
  const Matter & fluid = * mom->fluid();

  assert(mu_t->domain() == uvw.domain() );

  for_vijk( (*mu_t), i, j, k) {

    real ss;
    real ssd;

    /*----------------------+ 
    |  compute ss  [1/s^2]  |
    +----------------------*/
     {
      /* shear [1/s] */
      real s11, s12, s13, s21, s22, s23, s31, s32, s33;

      uvw.shear(i,j,k, &s11, &s12, &s13, 
                       &s21, &s22, &s23, 
                       &s31, &s32, &s33);

      /* 1/s^2 */
      ss = s11*s11 + s12*s12 + s13*s13 +
           s21*s21 + s22*s22 + s23*s23 +
           s31*s31 + s32*s32 + s33*s33;  
     }
  
    /*--------------+ 
    |  compute ssd  |
    +--------------*/
     {
      /* derivatives [1/s] */
      real du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
      uvw.grad(i,j,k, &du_dx, &du_dy, &du_dz,
                      &dv_dx, &dv_dy, &dv_dz,
                      &dw_dx, &dw_dy, &dw_dz);

      /* compute g_ij^2 [1/s^2] */
      real g11 = du_dx*du_dx + du_dy*dv_dx + du_dz*dw_dx;
      real g12 = du_dx*du_dy + du_dy*dv_dy + du_dz*dw_dy;
      real g13 = du_dx*du_dz + du_dy*dv_dz + du_dz*dw_dz;

      real g21 = dv_dx*du_dx + dv_dy*dv_dx + dv_dz*dw_dx;
      real g22 = dv_dx*du_dy + dv_dy*dv_dy + dv_dz*dw_dy;
      real g23 = dv_dx*du_dz + dv_dy*dv_dz + dv_dz*dw_dz;

      real g31 = dw_dx*du_dx + dw_dy*dv_dx + dw_dz*dw_dx;
      real g32 = dw_dx*du_dy + dw_dy*dv_dy + dw_dz*dw_dy;
      real g33 = dw_dx*du_dz + dw_dy*dv_dz + dw_dz*dw_dz;

      real s11d = 0.5*(g11 + g11) - 1.0/3.0 * g11;
      real s12d = 0.5*(g12 + g21);
      real s13d = 0.5*(g13 + g31);

      real s21d = 0.5*(g21 + g12);
      real s22d = 0.5*(g22 + g22) - 1.0/3.0 * g22;
      real s23d = 0.5*(g23 + g32);

      real s31d = 0.5*(g31 + g13);
      real s32d = 0.5*(g32 + g23);
      real s33d = 0.5*(g33 + g33) - 1.0/3.0 * g33;

      /* [1/s^4] */
      ssd = s11d*s11d + s12d*s12d + s13d*s13d +
            s21d*s21d + s22d*s22d + s23d*s23d +
            s31d*s31d + s32d*s32d + s33d*s33d;  
    }

    /*------------------------------------+
    |  delta = (dx * dy * dz)^(1/3)  [m]  |
    +------------------------------------*/
    real delta = pow((mu_t->dxc(i)*mu_t->dyc(j)*mu_t->dzc(k)),1.0/3.0); 

    /*------------------------------+
    |  compute turbulent viscosity  |
    +------------------------------*/
    real ssd_3_over_2 =       sqrt( ssd*ssd*ssd );           /* [1/s^6] */
    real ssd_5_over_4 = sqrt( sqrt( ssd*ssd*ssd*ssd*ssd ) ); /* [1/s^5] */ 
    real ss__5_over_2 =       sqrt( ss*ss*ss*ss*ss );        /* [1/s^5] */
    
    /* [kg/(ms)] */
    (*mu_t)[i][j][k] = fluid.rho(i,j,k) 
                     * c_w * c_w * delta * delta
                     * ssd_3_over_2 / (ss__5_over_2+ssd_5_over_4+boil::pico); 
  }

  /*----------------------------------+ 
  |  exchange buffers. important to   |
  |  exchange in the corners as well  |
  +----------------------------------*/
  mu_t->exchange_all();
}
