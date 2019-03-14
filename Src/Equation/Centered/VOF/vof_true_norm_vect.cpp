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

    mx[i][j][k] = nnx;
    my[i][j][k] = nny;
    mz[i][j][k] = nnz;
  }

  //boil::plot->plot(mx,my,mz, "mx-my-mz", time->current_step());
  //exit(0);

  mx.exchange_all();
  my.exchange_all();
  mz.exchange_all();

  return;
}
