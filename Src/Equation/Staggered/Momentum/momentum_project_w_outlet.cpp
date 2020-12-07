#include "momentum.h"

/******************************************************************************/
void Momentum::project_w_outlet(const Scalar & frc) {
  project_w_outlet(frc,u);

  return;
}

/******************************************************************************/
void Momentum::project_w_outlet(const Scalar & frc, Vector & veloc) {

  project(frc,u);

  real rho;
  real dt = time->dt();
  /* boundary */
  for_m(m)
    for(int b=0; b<veloc.bc(m).count(); b++) {
      if(veloc.bc(m).exists(b)) {
        if(veloc.bc(m).type(b) == BndType::outlet()) {
          Dir d = veloc.bc(m).direction(b);

          /* u */
          if(m == Comp::u() && (d == Dir::imin()||d == Dir::imax())) {
            for_vijk(veloc.bc(m).at(b),i,j,k) {
              rho = fluid()->rho(m,i,j,k);
              veloc[m][i][j][k] += (frc[i-1][j][k]-frc[i][j][k])/rho*dt
                                   /veloc.dxc(m,i);
            }
          }

          /* v */
          if(m == Comp::v() && (d == Dir::jmin()||d == Dir::jmax())) {
            for_vijk(veloc.bc(m).at(b),i,j,k) {
              rho = fluid()->rho(m,i,j,k);
              veloc[m][i][j][k] += (frc[i][j-1][k]-frc[i][j][k])/rho*dt
                                   /veloc.dyc(m,j);
            }
          }

          /* w */
          if(m == Comp::w() && (d == Dir::kmin()||d == Dir::kmax())) {
            for_vijk(veloc.bc(m).at(b),i,j,k) {
              rho = fluid()->rho(m,i,j,k);
              veloc[m][i][j][k] += (frc[i][j][k-1]-frc[i][j][k])/rho*dt
                                   /veloc.dzc(m,k);
            }
          }

        } /* if outlet */
      } /* if boundary exists */
    } /* m & b */
  
  veloc.exchange_all();

  return;
}
