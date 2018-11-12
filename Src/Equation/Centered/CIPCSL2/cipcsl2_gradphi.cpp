#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::gradphi(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector of grad(csa) at node point.
*         Resluts: nx, ny, nz
*******************************************************************************/

  real ni,nj,nk,magn;

  /* node base */
  for(int i=sca.si(); i<=sca.ei()+1; i++) {
    for(int j=sca.sj(); j<=sca.ej()+1; j++) {
      for(int k=sca.sk(); k<=sca.ek()+1; k++) {
        real dx=dxw(i);
        real dy=dys(j);
        real dz=dzb(k);
        nx[i][j][k] = 0.25*( (sca[i][j  ][k  ]-sca[i-1][j  ][k  ])/dx
                            +(sca[i][j-1][k  ]-sca[i-1][j-1][k  ])/dx
                            +(sca[i][j  ][k-1]-sca[i-1][j  ][k-1])/dx
                            +(sca[i][j-1][k-1]-sca[i-1][j-1][k-1])/dx);
        ny[i][j][k] = 0.25*( (sca[i  ][j][k  ]-sca[i  ][j-1][k  ])/dy
                            +(sca[i-1][j][k  ]-sca[i-1][j-1][k  ])/dy
                            +(sca[i  ][j][k-1]-sca[i  ][j-1][k-1])/dy
                            +(sca[i-1][j][k-1]-sca[i-1][j-1][k-1])/dy);
        nz[i][j][k] = 0.25*( (sca[i  ][j  ][k]-sca[i  ][j  ][k-1])/dz
                            +(sca[i  ][j-1][k]-sca[i  ][j-1][k-1])/dz
                            +(sca[i-1][j  ][k]-sca[i-1][j  ][k-1])/dz
                            +(sca[i-1][j-1][k]-sca[i-1][j-1][k-1])/dz);
      }
    }
  }
  
  insert_bc_gradphi(sca);

  for(int i=sca.si(); i<=sca.ei()+1; i++) {
    for(int j=sca.sj(); j<=sca.ej()+1; j++) {
      for(int k=sca.sk(); k<=sca.ek()+1; k++) {
        ni = nx[i][j][k];
        nj = ny[i][j][k];
        nk = nz[i][j][k];
        normalize(ni,nj,nk);
        nx[i][j][k] = ni;
        ny[i][j][k] = nj;
        nz[i][j][k] = nk;
      }
    }
  }

  return;
}
