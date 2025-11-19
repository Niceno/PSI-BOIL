#include "vof.h"

/******************************************************************************/
void VOF::standardized_norm_vect() {
/***************************************************************************//**
*  \brief Calculate normal vector at cell center.
*         Results: nx, ny, nz -- standardized normal vector 
*******************************************************************************/

  /* cell centered base, second order */
  for_aijk(i,j,k) {
    real mmx = mx[i][j][k];
    real mmy = my[i][j][k];
    real mmz = mz[i][j][k];

    real dnx = phi.dxc(i);
    real dny = phi.dyc(j);
    real dnz = phi.dzc(k);
    real nnx, nny, nnz;
    if(dnx==0.0||dny==0.0||dnz==0.0) {
      nnx = mmx;
      nny = mmy;
      nnz = mmz;
    } else {
      nnx = mmx*phi.dxc(i);
      nny = mmy*phi.dyc(j);
      nnz = mmz*phi.dzc(k);
    }

    normalize(nnx,nny,nnz);
 
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
