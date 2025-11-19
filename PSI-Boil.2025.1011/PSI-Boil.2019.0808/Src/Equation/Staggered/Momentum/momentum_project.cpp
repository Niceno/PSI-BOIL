#include "momentum.h"

/******************************************************************************/
void Momentum::project(const Scalar & frc) {
  project(frc,u);

  return;
}

/******************************************************************************/
void Momentum::project(const Scalar & frc, Vector & veloc) {
 
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
    if(ifull)
      for_mijk(m,i,j,k) {
        rho = fluid()->rho(m,i,j,k);
        veloc[m][i][j][k] += (frc[i-1][j][k]-frc[i][j][k]) / rho * time->dt()
                       / (veloc.dxc(m,i));
      }
    m = Comp::v();
    if(jfull)
      for_mijk(m,i,j,k) {
        rho = fluid()->rho(m,i,j,k);
        veloc[m][i][j][k] += (frc[i][j-1][k]-frc[i][j][k]) / rho * time->dt()
                       / (veloc.dyc(m,j));
      }
    m = Comp::w();
    if(kfull)
      for_mijk(m,i,j,k) {
        rho = fluid()->rho(m,i,j,k);
        veloc[m][i][j][k] += (frc[i][j][k-1]-frc[i][j][k]) / rho * time->dt()
                       / (veloc.dzc(m,k));
      }
  } else {
    m = Comp::u();
    if(ifull)
      for_vmijk(veloc,m,i,j,k) {
        if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i-1,j,k) )
          veloc[m][i][j][k] += (frc[i-1][j][k]-frc[i][j][k]) * time->dt() 
                          / fluid()->rho(m,i,j,k)
                          / (veloc.dxc(m,i));
      }
    m = Comp::v();
    if(jfull)
      for_vmijk(veloc,m,i,j,k) {
        if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i,j-1,k) )
          veloc[m][i][j][k] += (frc[i][j-1][k]-frc[i][j][k]) * time->dt()
                          / fluid()->rho(m,i,j,k)
                          / (veloc.dyc(m,j));
      }
    m = Comp::w();
    if(kfull)
      for_vmijk(veloc,m,i,j,k) {
        if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i,j,k-1) )
          veloc[m][i][j][k] += (frc[i][j][k-1]-frc[i][j][k]) * time->dt()
                          / fluid()->rho(m,i,j,k)
                          / (veloc.dzc(m,k));
      }
  }

  veloc.bnd_update_nooutlet();
  veloc.exchange_all();

  return;
}
