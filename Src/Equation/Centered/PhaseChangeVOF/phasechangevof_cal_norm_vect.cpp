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

    real nnx = mmx/phi.dxc(i);
    real nny = mmy/phi.dyc(j);
    real nnz = mmz/phi.dzc(k);

    real nnorm = nnx*nnx+nny*nny+nnz*nnz;
    nnorm = sqrt(nnorm)+boil::pico;

    nnx /= nnorm;
    nny /= nnorm;
    nnz /= nnorm;
  
    nx[i][j][k] = nnx;
    ny[i][j][k] = nny;
    nz[i][j][k] = nnz;
  }

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  return;
}
