#include "vector.h"

/******************************************************************************/
real Vector::shear_magn(const int i, const int j, const int k) const {
/*-------------------------------------------------------------------+
|  computes the magnitude of the shear at specified scalar location  |
|--------------------------------------------------------------------|
|  shear = sqrt( 2 * s_ij * s_ij )                                   |
|  where:                                                            |
|  s_ij = 1/2 ( du_i/dx_j + du_j/dx_i )                              |
|-------------------------------------------------------------------*/

  /*------------+
  |  get shear  |
  +------------*/
  real s11, s12, s13, s21, s22, s23, s31, s32, s33;

  shear(i, j, k, &s11, &s12, &s13,
                 &s21, &s22, &s23,
                 &s31, &s32, &s33);
  
  /*---------------------+
  |  magnitude of shear  |
  +---------------------*/
  real s_abs = s11*s11 + s12*s12 + s13*s13 +
               s21*s21 + s22*s22 + s23*s23 +
               s31*s31 + s32*s32 + s33*s33;  

  s_abs = sqrt( 2.0 * s_abs );

  return s_abs;
}
