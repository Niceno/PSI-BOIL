#include "vector.h"

/******************************************************************************/
void Vector::grad(const int i, const int j, const int k, 
                  real * du_dx, real * du_dy, real * du_dz,
                  real * dv_dx, real * dv_dy, real * dv_dz,
                  real * dw_dx, real * dw_dy, real * dw_dz) const {

  const Comp u=Comp::u();
  const Comp v=Comp::v();
  const Comp w=Comp::w();

  * du_dx = (vec[u][i+1][j][k] - vec[u][i][j][k]) / dxc(i);
  * dv_dy = (vec[v][i][j+1][k] - vec[v][i][j][k]) / dyc(j);
  * dw_dz = (vec[w][i][j][k+1] - vec[w][i][j][k]) / dzc(k);

  * dv_dx = 0.5 * (   (vec[v][i+1][j]  [k] - vec[v][i-1][j]  [k]) 
                    + (vec[v][i+1][j+1][k] - vec[v][i-1][j+1][k]) )
                / ( dxw(i) + dxe(i) );
  * dw_dx = 0.5 * (   (vec[w][i+1][j][k]   - vec[w][i-1][j][k]) 
                    + (vec[w][i+1][j][k+1] - vec[w][i-1][j][k+1]) )
                / ( dxw(i) + dxe(i) );

  * du_dy = 0.5 * (   (vec[u][i]  [j+1][k] - vec[u][i]  [j-1][k]) 
                    + (vec[u][i+1][j+1][k] - vec[u][i+1][j-1][k]) )
                / ( dys(j) + dyn(j) );
  * dw_dy = 0.5 * (   (vec[w][i][j+1][k]   - vec[w][i][j-1][k]  ) 
                    + (vec[w][i][j+1][k+1] - vec[w][i][j-1][k+1]) )
                / ( dys(j) + dyn(j) );

  * du_dz = 0.5 * (   (vec[u][i]  [j][k+1] - vec[u][i]  [j][k-1]) 
                    + (vec[u][i+1][j][k+1] - vec[u][i+1][j][k-1]) )
                / ( dzb(k) + dzt(k) );
  * dv_dz = 0.5 * (   (vec[v][i][j]  [k+1] - vec[v][i][j]  [k-1]) 
                    + (vec[v][i][j+1][k+1] - vec[v][i][j+1][k-1]) )
                / ( dzb(k) + dzt(k) );
}
