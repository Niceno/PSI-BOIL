#include "vector.h"

/***************************************************************************//**
*  Computes shear of velocity at the specified scalar location. Unit is [1/s] 
*******************************************************************************/
void Vector::shear(const int i, const int j, const int k, 
                   real * s11, real * s12, real * s13,
                   real * s21, real * s22, real * s23,
                   real * s31, real * s32, real * s33) const {
/*----------------------------------------------------------+
|  computes shear of velocity at specified scalar location  |
|-----------------------------------------------------------|
|  s_ij = 1/2 ( du_i/dx_j + du_j/dx_i )                     |
|----------------------------------------------------------*/

  real du_dx, du_dy, du_dz,
       dv_dx, dv_dy, dv_dz,
       dw_dx, dw_dy, dw_dz;

  /*-----------------------------+ 
  |  compute velocity gradients  |
  +-----------------------------*/
  grad( i, j, k, &du_dx, &du_dy, &du_dz,
                 &dv_dx, &dv_dy, &dv_dz,
                 &dw_dx, &dw_dy, &dw_dz);
 
  /*----------------+
  |  compute shear  |
  +----------------*/
  *s11 = 0.5 * (du_dx + du_dx);
  *s12 = 0.5 * (du_dy + dv_dx);
  *s13 = 0.5 * (du_dz + dw_dx);
 
  *s21 = *s12;
  *s22 = 0.5 * (dv_dy + dv_dy);
  *s23 = 0.5 * (dv_dz + dw_dy);
  
  *s31 = *s13;
  *s32 = *s23;
  *s33 = 0.5 * (dw_dz + dw_dz);
  
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_shear.cpp,v 1.3 2011/05/27 11:43:48 niceno Exp $'/
+-----------------------------------------------------------------------------*/
