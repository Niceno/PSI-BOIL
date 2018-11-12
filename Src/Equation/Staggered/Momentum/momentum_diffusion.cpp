#include "momentum.h"

/******************************************************************************/
void Momentum::diffusion() {
/*------------+ 
|       1  n  |
|  f += - D   |
|       2     |
+------------*/

  /* get time stepping coefficient */
  real tscn = diff_ts.N();
  real tsco = diff_ts.Nm1();
  assert( tscn > 0.0 );
  real tsc = tsco/tscn; /* 0.5 for c.n. 0.0 for fully implicit */

  for_m(m) {

    /*------------------------+ 
    |  x direction (w and e)  |
    +------------------------*/
    for_mijk(m,i,j,k) 
      fold[m][i][j][k] += (tsc * A[~m]->w[i][j][k] * (u[m][i-1][j][k] - u[m][i][j][k]) 
                         + tsc * A[~m]->e[i][j][k] * (u[m][i+1][j][k] - u[m][i][j][k])); 
    
    /*------------------------+ 
    |  y direction (s and n)  |
    +------------------------*/
    for_mijk(m,i,j,k) 
      fold[m][i][j][k] += (tsc * A[~m]->s[i][j][k] * (u[m][i][j-1][k] - u[m][i][j][k]) 
                         + tsc * A[~m]->n[i][j][k] * (u[m][i][j+1][k] - u[m][i][j][k])); 

    /*------------------------+ 
    |  z direction (b and t)  |
    +------------------------*/
    for_mijk(m,i,j,k) 
      fold[m][i][j][k] += (tsc * A[~m]->b[i][j][k] * (u[m][i][j][k-1] - u[m][i][j][k])
                         + tsc * A[~m]->t[i][j][k] * (u[m][i][j][k+1] - u[m][i][j][k]));
  }
}
