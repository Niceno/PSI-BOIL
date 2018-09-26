#include "levelset.h"

/******************************************************************************/
void LevelSet::gradphic() {
/***************************************************************************//**
*  \brief Calculate grad(csa)/|grad(csa)| at cell center.
*         Resluts: nx, ny, nz
*******************************************************************************/

  real ni,nj,nk,magn;

  /* cell centered base */
  for_ijk(i,j,k) {
    if(fabs(dflag[i][j][k])>real(nlayer-4)){
      nx[i][j][k] = 0.0;
      ny[i][j][k] = 0.0;
      nz[i][j][k] = 0.0;
    } else {
      nx[i][j][k] = (phi[i+1][j][k]-phi[i-1][j][k])/(dxw(i)+dxe(i));
      ny[i][j][k] = (phi[i][j+1][k]-phi[i][j-1][k])/(dys(j)+dyn(j));
      nz[i][j][k] = (phi[i][j][k+1]-phi[i][j][k-1])/(dzb(k)+dzt(k));
    }
  }

  insert_bc_gradphic(phi); 
  insert_bc_norm();

  for_avijk(phi,i,j,k) {
    if(fabs(dflag[i][j][k])>real(nlayer-4)){
      nx[i][j][k] = 0.0;
      ny[i][j][k] = 0.0;
      nz[i][j][k] = 0.0;
    } else {
      ni = nx[i][j][k];
      nj = ny[i][j][k];
      nk = nz[i][j][k];
      real magn = sqrt(ni*ni + nj*nj + nk*nk) + epsnorm;
      nx[i][j][k] = ni/magn;
      ny[i][j][k] = nj/magn;
      nz[i][j][k] = nk/magn;
    }
  }

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //boil::plot->plot(phi, "dist", time->current_step());
  //exit(0);

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: levelset_gradphic.cpp,v 1.2 2012/09/13 08:42:26 niceno Exp $'/
+-----------------------------------------------------------------------------*/
