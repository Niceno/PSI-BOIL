#include "levelset.h"
#include <iomanip>
//#define MARCHING

/******************************************************************************/
void LevelSet::tension(Vector * vec, const Matter matt, Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate surface tension
*         Algorithm
*           1st step: calculate curvature from distance function
*           2nd step: calculate body force
*         Variables
*           distance function : phi
*           curvature         : kappa
*           body force        : vec
*           color function    : sca
*******************************************************************************/

  boil::timer.start("levelsetdit tension");

  /*-----------+
  |  2nd step  |
  +-----------*/
  /* calculate curvature */
  curv();
  //boil::plot->plot(nx,ny,nz, "nx-ny-nz", time->current_step());

  /* wall adhesion */
  theta = cangle/ 180.0 * pi;
  bdcurv(phi,theta);
  insert_bc_kappa2(kappa);
  kappa.exchange();
  //boil::plot->plot(clr,kappa,dist, "clr-kappa-dist", time->current_step());
  //exit(0);

  /*-----------+
  |  3rd step  |
  +-----------*/
  real rho_diff = fabs(matt.rho(1)-matt.rho(0));
  real rho_ave = 0.5*(matt.rho(1)+matt.rho(0));
#ifdef MARCHING
  if(rho_diff!=0.0){
    gradphic(phi);
  }
#endif

  /* calculate body force */
  Comp m;
  m = Comp::u();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i-1][j][k]+kappa[i][j][k])
                         * (sca[i][j][k]-sca[i-1][j][k])/vec->dxc(m,i)
                         * vec->dV(m,i,j,k);
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i-1][j][k]+kappa[i][j][k])
#ifdef MARCHING
                         * marching_cube(m,i,j,k) / vec->dV(m,i,j,k)
                         * 0.5*(nx[i-1][j][k]+nx[i][j][k])
#else
                         * (matt.rho(i,j,k)-matt.rho(i-1,j,k))/vec->dxc(m,i)
                         / rho_diff * matt.rho(m,i,j,k) / rho_ave
#endif
                       * vec->dV(m,i,j,k);
    }
  }

  m = Comp::v();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j-1][k]+kappa[i][j][k])
                         * (sca[i][j][k] - sca[i][j-1][k])/vec->dxc(m,i)
                         * vec->dV(m,i,j,k);
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j-1][k]+kappa[i][j][k])
#ifdef MARCHING
                         * marching_cube(m,i,j,k) / vec->dV(m,i,j,k)
                         * 0.5*(ny[i][j-1][k]+ny[i][j][k])
#else
                         * (matt.rho(i,j,k)-matt.rho(i,j-1,k))/vec->dyc(m,j)
                         / rho_diff * matt.rho(m,i,j,k) / rho_ave
#endif
                         * vec->dV(m,i,j,k);
    }
  }

  m = Comp::w();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j][k-1]+kappa[i][j][k])
                         * (sca[i][j][k] - sca[i][j][k-1])/vec->dxc(m,i)
                         * vec->dV(m,i,j,k);
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j][k-1]+kappa[i][j][k])
#ifdef MARCHING
                         * marching_cube(m,i,j,k) / vec->dV(m,i,j,k)
                         * 0.5*(nz[i][j][k-1]+nz[i][j][k])
#else
                         * (matt.rho(i,j,k)-matt.rho(i,j,k-1))/vec->dzc(m,k)
                         / rho_diff * matt.rho(m,i,j,k) / rho_ave
#endif
                         * vec->dV(m,i,j,k);
    }
  }
  vec->exchange();

  boil::timer.stop("levelsetdit tension");
}
