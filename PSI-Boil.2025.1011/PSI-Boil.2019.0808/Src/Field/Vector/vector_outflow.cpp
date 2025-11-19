#include "vector.h"

/******************************************************************************/
real Vector::outflow(const int i, const int j, const int k) const {
/*-----------------------------------------------------------------------+
|  calculate volumetric outflow from a scalar cell: integral divergence  |
+-----------------------------------------------------------------------*/

  return - domain()->dSx(Sign::neg(),i,j,k)*vec[Comp::u()][i]  [j]  [k]
         + domain()->dSx(Sign::pos(),i,j,k)*vec[Comp::u()][i+1][j]  [k]
         - domain()->dSy(Sign::neg(),i,j,k)*vec[Comp::v()][i]  [j]  [k]
         + domain()->dSy(Sign::pos(),i,j,k)*vec[Comp::v()][i]  [j+1][k]
         - domain()->dSz(Sign::neg(),i,j,k)*vec[Comp::w()][i]  [j]  [k]
         + domain()->dSz(Sign::pos(),i,j,k)*vec[Comp::w()][i]  [j]  [k+1];
}
