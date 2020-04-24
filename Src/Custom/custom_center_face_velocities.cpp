#include "custom.h"
  
namespace boil {
  /******************************************************************************/
  void cell_center_velocities(const Vector & uvw,
                              Scalar & u, Scalar & v, Scalar & w) {
  /***************************************************************************//**
   \brief Interpolate velocity to cell center.
  *******************************************************************************/
    for_vijk(u,i,j,k) {
      u[i][j][k] = 0.5*(uvw[Comp::u()][i][j][k]+uvw[Comp::u()][i+1][j][k]);
      v[i][j][k] = 0.5*(uvw[Comp::v()][i][j][k]+uvw[Comp::v()][i][j+1][k]);
      w[i][j][k] = 0.5*(uvw[Comp::w()][i][j][k]+uvw[Comp::w()][i][j][k+1]);
    }

    return;
  }

  /******************************************************************************/
  void staggered_velocities(const Scalar & u, const Scalar & v, const Scalar & w,
                            Vector & uvw) {
  /***************************************************************************//**
   \brief Interpolate velocity to cell face (boundaries excluded).
  *******************************************************************************/
    Comp m = Comp::u();
    for_vmijk(uvw,m,i,j,k) {
      uvw[m][i][j][k] = 0.5*(u[i][j][k]+u[i-1][j][k]);
    }
    m = Comp::v();
    for_vmijk(uvw,m,i,j,k) {
      uvw[m][i][j][k] = 0.5*(v[i][j][k]+v[i][j-1][k]);
    }
    m = Comp::w();
    for_vmijk(uvw,m,i,j,k) {
      uvw[m][i][j][k] = 0.5*(w[i][j][k]+w[i][j][k-1]);
    }

    return;
  }
}
