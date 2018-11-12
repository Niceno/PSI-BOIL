#include "colorcip.h"

/******************************************************************************/
void ColorCIP::gradphi(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate grad(csa)/|grad(csa)| at node point.
*         Reference: J.U.Brackbill, et al.,"A Continum method for modeling
*                    surface tension",J.Comp.phys.,Vol.100,pp.335-354,1992
*                    Equation (41)
*         Resluts: nx, ny, nz
*******************************************************************************/

  real ni,nj,nk,magn;

  /* node base */
  for(int i=1; i<=sca.ei()+1; i++) {
    for(int j=1; j<=sca.ej()+1; j++) {
      for(int k=1; k<=sca.ek()+1; k++) {
        real dx=dxw(i);
        real dy=dys(j);
        real dz=dzb(k);
        ni = 0.25*( (sca[i][j  ][k  ]-sca[i-1][j  ][k  ])/dx
                   +(sca[i][j-1][k  ]-sca[i-1][j-1][k  ])/dx
                   +(sca[i][j  ][k-1]-sca[i-1][j  ][k-1])/dx
                   +(sca[i][j-1][k-1]-sca[i-1][j-1][k-1])/dx);
        nj = 0.25*( (sca[i  ][j][k  ]-sca[i  ][j-1][k  ])/dy
                   +(sca[i-1][j][k  ]-sca[i-1][j-1][k  ])/dy
                   +(sca[i  ][j][k-1]-sca[i  ][j-1][k-1])/dy
                   +(sca[i-1][j][k-1]-sca[i-1][j-1][k-1])/dy);
        nk = 0.25*( (sca[i  ][j  ][k]-sca[i  ][j  ][k-1])/dz
                   +(sca[i  ][j-1][k]-sca[i  ][j-1][k-1])/dz
                   +(sca[i-1][j  ][k]-sca[i-1][j  ][k-1])/dz
                   +(sca[i-1][j-1][k]-sca[i-1][j-1][k-1])/dz);
        real magn = sqrt(ni*ni + nj*nj + nk*nk) + epsnorm;

        nx[i][j][k] = ni;
        ny[i][j][k] = nj;
        nz[i][j][k] = nk;
        nmag[i][j][k] = magn;
      }
    }
  }
  return;
}
