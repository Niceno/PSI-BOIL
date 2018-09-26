#include "phasechange.h"

/******************************************************************************/
void PhaseChange::gradphic(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate grad(csa)/|grad(csa)| at cell center.
*         Resluts: nx, ny, nz
*******************************************************************************/

  /* cell centered base */
  for_ijk(i,j,k) {
    if(abs(dflag[i][j][k])>real(nlayer-4)){
      nx[i][j][k] = 0.0;
      ny[i][j][k] = 0.0;
      nz[i][j][k] = 0.0;
    } else {
      nx[i][j][k] = (sca[i+1][j][k]-sca[i-1][j][k])/(dxw(i)+dxe(i));
      ny[i][j][k] = (sca[i][j+1][k]-sca[i][j-1][k])/(dys(j)+dyn(j));
      nz[i][j][k] = (sca[i][j][k+1]-sca[i][j][k-1])/(dzb(k)+dzt(k));
    }
  }

  insert_bc_gradphic(sca); 
  insert_bc_norm();

  /* normalize */
  for_avijk(sca,i,j,k) {
    if(abs(dflag[i][j][k])>real(nlayer-4)){
      nx[i][j][k] = 0.0;
      ny[i][j][k] = 0.0;
      nz[i][j][k] = 0.0;
    } else {
      normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }
  }

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //boil::plot->plot(sca, "dist", time->current_step());
  //exit(0);

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: phasechange_gradphic.cpp,v 1.3 2015/05/05 14:50:51 sato Exp $'/
+-----------------------------------------------------------------------------*/
