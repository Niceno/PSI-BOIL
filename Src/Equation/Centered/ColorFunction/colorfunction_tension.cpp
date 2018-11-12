#include "colorfunction.h"

/******************************************************************************/
void ColorFunction::tension(Vector * vec, const Matter matt) {

  /* limit surface tension on one part of the domain only: */
  /* 0.01 -> 0.35 proved good for PressureJump */
  for_ijk(i,j,k) {
    if( phi[i][j][k] > 0.01 && phi[i][j][k] < 0.50 ) {
      surf[i][j][k] = 1.0;
    } else {
      surf[i][j][k] = 0.0;
    }
  }
  for_ijk(i,j,k) {
    surf[i][j][k] = 4.0*phi[i][j][k]*(1.0 - phi[i][j][k]);
  }
  surf.exchange();

  Comp m;

  m = Comp::u();
  for_vmijk((*vec),m,i,j,k) {
    real nw = nx[i-1][j][k];
    real ne = nx[i]  [j][k];
    /*.................... i ............... i-1 ....................*/
    real ns = 0.25 * (  ny[i][j]  [k]   + ny[i-1][j]  [k]     /* j   */
                      + ny[i][j-1][k]   + ny[i-1][j-1][k]  ); /* j-1 */
    real nn = 0.25 * (  ny[i][j]  [k]   + ny[i-1][j]  [k]     /* j   */
                      + ny[i][j+1][k]   + ny[i-1][j+1][k]  ); /* j+1 */
    real nb = 0.25 * (  nz[i][j]  [k]   + nz[i-1][j]  [k]     /* k   */
                      + nz[i][j]  [k-1] + nz[i-1][j]  [k-1]); /* k-1 */
    real nt = 0.25 * (  nz[i][j]  [k]   + nz[i-1][j]  [k]     /* k   */
                      + nz[i][j]  [k+1] + nz[i-1][j]  [k+1]); /* k+1 */

    (*vec)[m][i][j][k] = - 0.5 * (surf[i-1][j][k] + surf[i][j][k]) // on surface?
                       * (phi[i][j][k] - phi[i-1][j][k])/vec->dxc(m,i) // grad(phi)
                       * vec->dV(m,i,j,k)                              // dvol
                       * matt.sigma(m,i,j,k)                           // sigma
                       * (   (ne-nw)/vec->dxc(m,i)                     // nabla(n)
                           + (nn-ns)/vec->dyc(m,j) 
                           + (nt-nb)/vec->dzc(m,k) );
  }

  m = Comp::v();
  for_vmijk((*vec),m,i,j,k) {
    real ns = ny[i][j-1][k];
    real nn = ny[i][j]  [k];
    /*......................... j ............... j-1 ...............*/
    real nw = 0.25 * (  nx[i]  [j][k]   + nx[i]  [j-1][k]     /* i   */
                      + nx[i-1][j][k]   + nx[i-1][j-1][k]  ); /* i-1 */
    real ne = 0.25 * (  nx[i]  [j][k]   + nx[i]  [j-1][k]     /* i   */
                      + nx[i+1][j][k]   + nx[i+1][j-1][k]  ); /* i+1 */
    real nb = 0.25 * (  nz[i]  [j][k]   + nz[i]  [j-1][k]     /* k   */
                      + nz[i]  [j][k-1] + nz[i]  [j-1][k-1]); /* k-1 */
    real nt = 0.25 * (  nz[i]  [j][k]   + nz[i]  [j-1][k]     /* k   */
                      + nz[i]  [j][k+1] + nz[i]  [j-1][k+1]); /* k+1 */

    (*vec)[m][i][j][k] = - 0.5 * (surf[i][j-1][k] + surf[i][j][k]) // on surface?
                       * (phi[i][j][k] - phi[i][j-1][k])/vec->dyc(m,j) // grad(phi)
                       * vec->dV(m,i,j,k)                              // dvol
                       * matt.sigma(m,i,j,k)                           // sigma
                       * (   (ne-nw)/vec->dxc(m,i)                     // nabla(n)
                           + (nn-ns)/vec->dyc(m,j) 
                           + (nt-nb)/vec->dzc(m,k) );
  }

  m = Comp::w();
  for_vmijk((*vec),m,i,j,k) {
    /*.............................. k ..............  k-1 ..........*/
    real nw = 0.25 * (  nx[i]  [j]  [k] + nx[i]  [j]  [k-1]   /* i   */
                      + nx[i-1][j]  [k] + nx[i-1][j]  [k-1]); /* i-1 */
    real ne = 0.25 * (  nx[i]  [j]  [k] + nx[i]  [j]  [k-1]   /* i   */
                      + nx[i+1][j]  [k] + nx[i+1][j]  [k-1]); /* i+1 */
    real ns = 0.25 * (  ny[i]  [j]  [k] + ny[i]  [j]  [k-1]   /* j   */
                      + ny[i]  [j-1][k] + ny[i]  [j-1][k-1]); /* j-1 */
    real nn = 0.25 * (  ny[i]  [j]  [k] + ny[i]  [j]  [k-1]   /* j   */
                      + ny[i]  [j+1][k] + ny[i]  [j+1][k-1]); /* j+1 */
    real nb = nz[i][j][k-1];
    real nt = nz[i][j][k];

    (*vec)[m][i][j][k] = - 0.5 * (surf[i][j][k-1] + surf[i][j][k]) // on surface?
                       * (phi[i][j][k] - phi[i][j][k-1])/vec->dzc(m,k) // grad(phi)
                       * vec->dV(m,i,j,k)                              // dvol
                       * matt.sigma(m,i,j,k)                           // sigma
                       * (   (ne-nw)/vec->dxc(m,i)                     // nabla(n)
                           + (nn-ns)/vec->dyc(m,j) 
                           + (nt-nb)/vec->dzc(m,k) );
  }


  vec->exchange();
}
