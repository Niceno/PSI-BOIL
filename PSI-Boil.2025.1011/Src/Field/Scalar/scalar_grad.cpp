#include "scalar.h"

/******************************************************************************/
void Scalar::grad( const int i, const int j, const int k,
                   real * d_x, real * d_y, real * d_z) const {

  /* this is probably not the most accurate for immersed cells */
  * d_x = (val[i+1][j][k] - val[i-1][j][k]) / (dxe(i) + dxw(i));
  * d_y = (val[i][j+1][k] - val[i][j-1][k]) / (dys(j) + dyn(j));
  * d_z = (val[i][j][k+1] - val[i][j][k-1]) / (dzb(k) + dzt(k));
}

/******************************************************************************/
void Scalar::grad_abs( const int i, const int j, const int k,
                       real * d_x, real * d_y, real * d_z) const {

  /* this is probably not the most accurate for immersed cells */
  real dx = (val[i+1][j][k] - val[i-1][j][k]) / (dxe(i) + dxw(i));
  real dy = (val[i][j+1][k] - val[i][j-1][k]) / (dys(j) + dyn(j));
  real dz = (val[i][j][k+1] - val[i][j][k-1]) / (dzb(k) + dzt(k));

  real d = sqrt( dx*dx + dy*dy + dz*dz );

  *d_x = dx / d;  
  *d_y = dy / d;  
  *d_z = dz / d;  
}
