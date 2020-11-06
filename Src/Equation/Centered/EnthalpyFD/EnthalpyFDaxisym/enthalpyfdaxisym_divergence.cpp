#include "enthalpyfdaxisym.h"

/***************************************************************************//**
*  \brief Divergence of flux (to be overwritten by child class).
*******************************************************************************/
real EnthalpyFDaxisym::neg_div_x(const int i, const int j, const int k,
                                 const Vector & flux) {
  return (flux[Comp::u()][i  ][j][k]*flux.xc(Comp::u(),i  )
         -flux[Comp::u()][i+1][j][k]*flux.xc(Comp::u(),i+1))
         /flux.dxe(Comp::u(),i)/flux.xn(Comp::u(),i+1);
}

real EnthalpyFDaxisym::neg_div_y(const int i, const int j, const int k,
                                 const Vector & flux) {
  return 0.0;
}

real EnthalpyFDaxisym::neg_div_z(const int i, const int j, const int k,
                                 const Vector & flux) {
  return (flux[Comp::w()][i][j][k]-flux[Comp::w()][i][j][k+1])
         /flux.dzt(Comp::k(),k);
}
