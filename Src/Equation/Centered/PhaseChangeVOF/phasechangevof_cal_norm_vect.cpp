#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::cal_norm_vect() {
/***************************************************************************//**
*  \brief Calculate normal vector at cell center.
*         Results: nx, ny, nz
*******************************************************************************/

  /* cell centered base, second order */
  for_aijk(i,j,k) {
    real mmx = mx[i][j][k];
    real mmy = my[i][j][k];
    real mmz = mz[i][j][k];

#if 1
    real dnx = phi.dxc(i);
    real dny = phi.dyc(j);
    real dnz = phi.dzc(k);
    real nnx, nny, nnz;
    if(dnx==0.0||dny==0.0||dnz==0.0) {
      nnx = mmx;
      nny = mmy;
      nnz = mmz;
    } else {
      nnx = mmx/phi.dxc(i);
      nny = mmy/phi.dyc(j);
      nnz = mmz/phi.dzc(k);
    }
 
    real nnorm = nnx*nnx+nny*nny+nnz*nnz;
    nnorm = sqrt(nnorm)+boil::pico;

    nnx /= nnorm;
    nny /= nnorm;
    nnz /= nnorm;
  
#elif 0
    nnx = phi.xc(i)/sqrt(1.74505e-06*1.74505e-06-phi.xc(i)*phi.xc(i));
    nny = 0.0;
    nnz = -1.0;
    nnorm = nnx*nnx+nny*nny+nnz*nnz;
    nnorm = sqrt(nnorm)+boil::pico;
    nnx /= nnorm;
    nny /= nnorm;
    nnz /= nnorm;
#else
    real nnx = mmx;
    real nny = mmy;
    real nnz = mmz;
#endif

    nx[i][j][k] = nnx;
    ny[i][j][k] = nny;
    nz[i][j][k] = nnz;
  }

  //boil::plot->plot(nx,ny,nz, "nx-ny-nz", time->current_step());
  //exit(0);

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  return;
}
