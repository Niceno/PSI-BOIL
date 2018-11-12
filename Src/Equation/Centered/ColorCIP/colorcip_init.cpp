#include "colorcip.h"

/******************************************************************************/
void ColorCIP::init() {

  insert_bc(phi);
  phi.exchange_all();

  for_aijk(i,j,k)
    clr[i][j][k]=phi[i][j][k];

#ifdef USE_TAN
  for_aijk(i,j,k)
     clr[i][j][k]=tan((clr[i][j][k]-0.5)*tanfac*pi);
#endif

  for_aijk(i,j,k) {
    gpx[i][j][k]=0.0;
    gpy[i][j][k]=0.0;
    gpz[i][j][k]=0.0;
  }

  for_ijk(i,j,k) {
    gpx[i][j][k] = (clr[i+1][j][k]-clr[i-1][j][k])/(dxw(i)+dxe(i));
    gpy[i][j][k] = (clr[i][j+1][k]-clr[i][j-1][k])/(dys(j)+dyn(j));
    gpz[i][j][k] = (clr[i][j][k+1]-clr[i][j][k-1])/(dzb(k)+dzt(k));
  }

#ifdef USE_TAN
  for_aijk(i,j,k)
     clr[i][j][k]=atan(clr[i][j][k])/(tanfac*pi)+0.5;
#endif

  gpx.bnd_grad_update(Comp::i());
  gpy.bnd_grad_update(Comp::j());
  gpz.bnd_grad_update(Comp::k());

  gpx.exchange();
  gpy.exchange();
  gpz.exchange();

  return;
}
