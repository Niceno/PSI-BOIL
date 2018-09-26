#include "colorcip.h"

/******************************************************************************/
void ColorCIP::gradphic(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate grad(csa)/|grad(csa)| at cell center.
*         Resluts: nx, ny, nz
*******************************************************************************/

  real ni,nj,nk,magn;

  /* cell centered base */
  for_ijk(i,j,k) {
    ni = (sca[i+1][j][k]-sca[i-1][j][k])/(dxw(i)+dxe(i));
    nj = (sca[i][j+1][k]-sca[i][j-1][k])/(dys(j)+dyn(j));
    nk = (sca[i][j][k+1]-sca[i][j][k-1])/(dzb(k)+dzt(k));

    real magn = sqrt(ni*ni + nj*nj + nk*nk) + epsnorm;

    nx[i][j][k] = ni/magn;
    ny[i][j][k] = nj/magn;
    nz[i][j][k] = nk/magn;
  }

  nx.bnd_grad_update(Comp::i());
  ny.bnd_grad_update(Comp::j());
  nz.bnd_grad_update(Comp::k());
  nx.exchange();
  ny.exchange();
  nz.exchange();

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_gradphic.cpp,v 1.3 2014/10/15 13:39:30 niceno Exp $'/
+-----------------------------------------------------------------------------*/
