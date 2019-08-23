#include "vof.h"

/******************************************************************************/
void VOF::true_norm_vect() {
/***************************************************************************//**
*  \brief Calculate normal vector at cell center.
*         Results: mx, my, mz -- true normal vector 
*******************************************************************************/

  /* cell centered base, second order */
  for_aijk(i,j,k) {
    real mmx = nx[i][j][k];
    real mmy = ny[i][j][k];
    real mmz = nz[i][j][k];

    real dnx = nx.dxc(i);
    real dny = nx.dyc(j);
    real dnz = nx.dzc(k);
    real nnx, nny, nnz;
    if(dnx==0.0||dny==0.0||dnz==0.0) {
      nnx = mmx;
      nny = mmy;
      nnz = mmz;
    } else {
      nnx = mmx/nx.dxc(i);
      nny = mmy/nx.dyc(j);
      nnz = mmz/nx.dzc(k);
    }

    normalize(nnx,nny,nnz);
 
    mx[i][j][k] = nnx;
    my[i][j][k] = nny;
    mz[i][j][k] = nnz;

  }

  //boil::plot->plot(mx,my,mz, "mx-my-mz", time->current_step());
  //exit(0);

  mx.bnd_update();
  my.bnd_update();
  mz.bnd_update();

  mx.exchange_all();
  my.exchange_all();
  mz.exchange_all();

  return;
}
