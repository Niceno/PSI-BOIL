#include "momentum.h"

/******************************************************************************/
void Momentum::get_eps(Scalar * eps) {
/*-----------------------------------------------------------------------------+
|                                                                              |
|         |     du/dx         1/2(du/dy+dv/dx)  1/2(du/dz+dw/dx) |             |
|         |                                                      |             |
|  S_ij = | 1/2(du/dy+dv/dx)      dv/dy         1/2(dv/dz+dw/dy) |             |
|         |                                                      |             |
|         | 1/2(du/dz+dw/dx)  1/2(dv/dz+dw/dy)      dw/dz        |             |
|                                                                              |
|  eps = 2 * nu * S_ij * S_ij                                                  |
|                                                                              |
|  S_ij * S_ij = S_11^2 + S_22^2 + S_33^2                                      |
|              + S_12^2 + S_13^2                                               |
|              + S_21^2 + S_23^2                                               |
|              + S_31^2 + S_32^2                                               |
|                                                                              |
|  since S_ij = S_ji:                                                          |
|                                                                              |
|  S_ij * S_ij = S_11^2 + S_22^2 + S_33^2                                      |
|              + 2.0*S_12^2 + 2.0*S_13^2 + 2.0*S_23^2                          |
|                                                                              |
|  S_11^2 = (du/dx)^2                                                          |
|  S_22^2 = (dv/dy)^2                                                          |
|  S_33^2 = (dw/dz)^2                                                          |
|                                                                              |
|  S_12^2 = 1/4 * (du/dy+dv/dx)^2                                              |
|  S_13^2 = 1/4 * (du/dz+dw/dx)^2                                              |
|  S_23^2 = 1/4 * (dv/dz+dw/dy)^2                                              |
|                                                                              |
|  S_ij * S_ij = (du/dx)^2 + (dv/dy)^2 + (dw/dz)^2                             |
|              + 1/2 * (du/dy+dv/dx)^2                                         |
|              + 1/2 * (du/dz+dw/dx)^2                                         |
|              + 1/2 * (dv/dz+dw/dy)^2                                         |
|                                                                              |
|  eps = nu * { 2.0 * [ (du/dx)^2 + (dv/dy)^2 + (dw/dz)^2 ]                    |
|              + (du/dy+dv/dx)^2                                               |
|              + (du/dz+dw/dx)^2                                               |
|              + (dv/dz+dw/dy)^2 }                                             |
|                                                                              |
|  eps = nu * { 2.0 * [ (du/dx)^2 + (dv/dy)^2 + (dw/dz)^2 ]                    |
|              + du_dy2 + 2.0*(du/dy*dv/dx) + dv_dx2                           |
|              + du_dz2 + 2.0*(du/dz*dw/dx) + dw_dx2                           |
|              + dv_dz2 + 2.0*(dv/dz*dw/dy) + dw_dy2 }                         |
|                                                                              |
|  eps = nu * { 2.0 * [ (du/dx)^2 + (dv/dy)^2 + (dw/dz)^2 ]                    |
|              + du_dy2 + 2.0*(du/dy*dv/dx) + dv_dx2                           |
|              + du_dz2 + 2.0*(du/dz*dw/dx) + dw_dx2                           |
|              + dv_dz2 + 2.0*(dv/dz*dw/dy) + dw_dy2 }                         |
|                                                                              |
+-----------------------------------------------------------------------------*/

  u.exchange();
  
  for_vijk((*eps),i,j,k) {

    real du_dx = (u[Comp::u()][i+1][j][k] - u[Comp::u()][i][j][k]) / u.dxc(i);

    real du_dy = 0.25 
               * (   (u[Comp::u()][i]  [j+1][k] - u[Comp::u()][i]  [j-1][k]) 
                   + (u[Comp::u()][i+1][j+1][k] - u[Comp::u()][i+1][j-1][k]) )
               / u.dyc(j);

    real du_dz = 0.25 
               * (   (u[Comp::u()][i]  [j][k+1] - u[Comp::u()][i]  [j][k-1]) 
                   + (u[Comp::u()][i+1][j][k+1] - u[Comp::u()][i+1][j][k-1]) )
               / u.dzc(k);

    real dv_dx = 0.25 
               * (   (u[Comp::v()][i+1][j]  [k] - u[Comp::v()][i-1][j]  [k]) 
                   + (u[Comp::v()][i+1][j+1][k] - u[Comp::v()][i-1][j+1][k]) )
               / u.dxc(i);

    real dv_dy = (u[Comp::v()][i][j+1][k] - u[Comp::v()][i][j][k]) / u.dyc(j);

    real dv_dz = 0.25 
               * (   (u[Comp::v()][i][j]  [k+1] - u[Comp::v()][i][j]  [k-1]) 
                   + (u[Comp::v()][i][j+1][k+1] - u[Comp::v()][i][j+1][k-1]) )
               / u.dzc(k);

    real dw_dx = 0.25 
               * (   (u[Comp::w()][i+1][j][k]   - u[Comp::w()][i-1][j][k]) 
                   + (u[Comp::w()][i+1][j][k+1] - u[Comp::w()][i-1][j][k+1]) )
               / u.dxc(i);

    real dw_dy = 0.25 
               * (   (u[Comp::w()][i][j+1][k]   - u[Comp::w()][i][j-1][k]  ) 
                   + (u[Comp::w()][i][j+1][k+1] - u[Comp::w()][i][j-1][k+1]) )
               / u.dyc(j);

    real dw_dz = (u[Comp::w()][i][j][k+1] - u[Comp::w()][i][j][k]) / u.dzc(k);

    real du_dx2 = du_dx * du_dx;
    real du_dy2 = du_dy * du_dy;
    real du_dz2 = du_dz * du_dz;

    real dv_dx2 = dv_dx * dv_dx;
    real dv_dy2 = dv_dy * dv_dy;
    real dv_dz2 = dv_dz * dv_dz;

    real dw_dx2 = dw_dx * dw_dx;
    real dw_dy2 = dw_dy * dw_dy;
    real dw_dz2 = dw_dz * dw_dz;

    real mu = fluid()->mu(i,j,k);

    (*eps)[i][j][k] = mu * ( 2.0 * (du_dx2  + 
                                    dv_dy2  + 
                                    dw_dz2) +

                                    du_dy2 + du_dz2 +
                                    dv_dx2 + dv_dz2 +
                                    dw_dx2 + dw_dy2 +

                             2.0 * (du_dy * dv_dx + 
                                    du_dz * dw_dx + 
                                    dv_dz * dw_dy) );  
  }

  return;
}
