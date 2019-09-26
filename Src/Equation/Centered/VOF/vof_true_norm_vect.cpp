#include "vof.h"

/******************************************************************************/
void VOF::true_norm_vect() {
/***************************************************************************//**
*  \brief Calculate normal vector at cell center.
*         Results: mx, my, mz -- true normal vector 
*******************************************************************************/

  /* cell centered base, second order */
  for_aijk(i,j,k) {
    real nnx = nx[i][j][k];
    real nny = ny[i][j][k];
    real nnz = nz[i][j][k];

    real dnx = nx.dxc(i);
    real dny = nx.dyc(j);
    real dnz = nx.dzc(k);
    real mmx, mmy, mmz;
    if(dnx==0.0||dny==0.0||dnz==0.0) {
      mmx = nnx;
      mmy = nny;
      mmz = nnz;
    } else {
      mmx = nnx/nx.dxc(i);
      mmy = nny/nx.dyc(j);
      mmz = nnz/nx.dzc(k);
    }

    normalize(mmx,mmy,mmz);
 
    mx[i][j][k] = mmx;
    my[i][j][k] = mmy;
    mz[i][j][k] = mmz;

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
