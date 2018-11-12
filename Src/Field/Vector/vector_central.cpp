#include "vector.h"

/******************************************************************************/
void Vector::central( const int i, const int j, const int k,
                      real * u, real * v, real * w) const {

  /* this is not the most accurate for immersed cells 
     see how it is done in Plot object */
  * u = 0.5 * ( vec[Comp::u()][i+1][j][k] + vec[Comp::u()][i][j][k] );
  * v = 0.5 * ( vec[Comp::v()][i][j+1][k] + vec[Comp::v()][i][j][k] );
  * w = 0.5 * ( vec[Comp::w()][i][j][k+1] + vec[Comp::w()][i][j][k] );
}

/******************************************************************************/
void Vector::central( const int i, const int j, const int k,
                      real * uvw) const {

  /* this is not the most accurate for immersed cells 
     see how it is done in Plot object */
  uvw[0] = 0.5 * ( vec[Comp::u()][i+1][j][k] + vec[Comp::u()][i][j][k] );
  uvw[1] = 0.5 * ( vec[Comp::v()][i][j+1][k] + vec[Comp::v()][i][j][k] );
  uvw[2] = 0.5 * ( vec[Comp::w()][i][j][k+1] + vec[Comp::w()][i][j][k] );
}
