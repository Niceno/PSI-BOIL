#include "colorcip.h"
#include <iomanip>

/******************************************************************************/
void ColorCIP::smooth(const Scalar & sca, Scalar & scb, const int itnum) {
/*---------------------------------------------------+
| smooth and cut-off phi.                            |
| if itnum=0, then only the cut-off function works.  |
|   input:sca, itnum                                 |
|   output:scb                                       |
|   temporary:stmp                                   |
+---------------------------------------------------*/

  for_aijk(i,j,k)
    scb[i][j][k]=sca[i][j][k];

  if (itnum>=1) {
    /*------------+
    |   smooth    |
    +------------*/

    /* coefficients for smooth */
    const real c1 = 1.0/(6.0+12.0/sqrt(2.0)+8/sqrt(3.0));
    const real c2 = c1/sqrt(2.0);
    const real c3 = c1/sqrt(3.0);

    /* iterate */
    for(int it=0; it<itnum; it++) {
      for_ijk(i,j,k) {
        stmp[i][j][k] =0.5*scb[i][j][k] + 0.5/(1.0+6.0*c1+12.0*c2+8.0*c3)
                         *(scb[i][j][k]
                      +c1*(scb[i-1][j][k]+scb[i+1][j][k]
                          +scb[i][j-1][k]+scb[i][j+1][k]
                          +scb[i][j][k-1]+scb[i][j][k+1])
                      +c2*(scb[i-1][j-1][k]+scb[i-1][j+1][k] 
                          +scb[i+1][j-1][k]+scb[i+1][j+1][k]
                          +scb[i-1][j][k-1]+scb[i-1][j][k+1]
                          +scb[i+1][j][k-1]+scb[i+1][j][k+1]
                          +scb[i][j-1][k-1]+scb[i][j-1][k+1]
                          +scb[i][j+1][k-1]+scb[i][j+1][k+1])
                      +c3*(scb[i-1][j-1][k-1]+scb[i-1][j-1][k+1]
                          +scb[i-1][j+1][k-1]+scb[i-1][j+1][k+1]
                          +scb[i+1][j-1][k-1]+scb[i+1][j-1][k+1]
                          +scb[i+1][j+1][k-1]+scb[i+1][j+1][k+1]));
      }
      insert_bc(stmp);
      stmp.exchange_all();
      for_aijk(i,j,k)
        scb[i][j][k]=stmp[i][j][k];
    }
  }

  /*----------+
  |  cut-off  |
  +----------*/

  for_aijk(i,j,k)
    scb[i][j][k]=std::max(0.0,(std::min(1.0,scb[i][j][k])));

  return;
}
