#include "momentum.h"

/******************************************************************************/
void Momentum::project(const Scalar & frc) {
 
  Comp m;
  real rho;

  /*----------------------------------------------------+
  |  project the velocities on the computational cells  |
  +- - - - - - - - - - - - - - - - - - - - - - - - - - -+
  |     not the most efficient -> too many divisions    |
  +- - - - - - - - - - - - - - - - - - - - - - - - - - -+
  |      do not try to introduce some integral form     |
  |       here gradients are constant in each cell      |
  |          and depend on boundary values only         |
  +----------------------------------------------------*/

  if( dom->ibody().nccells() == 0 ) {
    m = Comp::u();
    for_mijk(m,i,j,k) {
      rho = fluid()->rho(m,i,j,k);
      u[m][i][j][k] += (frc[i-1][j][k]-frc[i][j][k]) / rho * time->dt()
                     / (u.dxc(m,i));
    }
    m = Comp::v();
    for_mijk(m,i,j,k) {
      rho = fluid()->rho(m,i,j,k);
      u[m][i][j][k] += (frc[i][j-1][k]-frc[i][j][k]) / rho * time->dt()
                     / (u.dyc(m,j));
    }
    m = Comp::w();
    for_mijk(m,i,j,k) {
      rho = fluid()->rho(m,i,j,k);
      u[m][i][j][k] += (frc[i][j][k-1]-frc[i][j][k]) / rho * time->dt()
                     / (u.dzc(m,k));
    }
  } else {
    m = Comp::u();
    for_vmijk(u,m,i,j,k) {
      if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i-1,j,k) )
        u[m][i][j][k] += (frc[i-1][j][k]-frc[i][j][k]) * time->dt() 
                        / fluid()->rho(m,i,j,k)
                        / (u.dxc(m,i));
    }
    m = Comp::v();
    for_vmijk(u,m,i,j,k) {
      if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i,j-1,k) )
        u[m][i][j][k] += (frc[i][j-1][k]-frc[i][j][k]) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dyc(m,j));
    }
    m = Comp::w();
    for_vmijk(u,m,i,j,k) {
      if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i,j,k-1) )
        u[m][i][j][k] += (frc[i][j][k-1]-frc[i][j][k]) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dzc(m,k));
    }
  }
  insert_bc();
  u.exchange_all();
}
