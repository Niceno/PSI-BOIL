#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Divergence of flux (to be overwritten by child class).
*******************************************************************************/
real EnthalpyFD::neg_div_x(const int i, const int j, const int k,
                           const Vector & flux) {
  return (flux[Comp::u()][i][j][k]-flux[Comp::u()][i+1][j][k])
         /flux.dxe(Comp::u(),i);
}

real EnthalpyFD::neg_div_y(const int i, const int j, const int k,
                           const Vector & flux) {
  return (flux[Comp::v()][i][j][k]-flux[Comp::v()][i][j+1][k])
         /flux.dyn(Comp::v(),j);
}

real EnthalpyFD::neg_div_z(const int i, const int j, const int k,
                           const Vector & flux) {
  return (flux[Comp::w()][i][j][k]-flux[Comp::w()][i][j][k+1])
         /flux.dzt(Comp::k(),k);
}
