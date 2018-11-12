#include "levelset.h"

/******************************************************************************/
void LevelSet::curv() {
/***************************************************************************//**
*  \brief Calculate curvature.
*     output: kappa
*******************************************************************************/

#if 0
  /* node base */
  gradphi();
  for_ijk(i,j,k) {
    real dx=clr.dxc(i);
    real dy=clr.dyc(j);
    real dz=clr.dzc(k);
    real nxcent = 0.125 * (nx[i  ][j  ][k  ]+nx[i+1][j  ][k  ]
                          +nx[i  ][j+1][k  ]+nx[i+1][j+1][k  ]
                          +nx[i  ][j  ][k+1]+nx[i+1][j  ][k+1]
                          +nx[i  ][j+1][k+1]+nx[i+1][j+1][k+1]);
    real nycent = 0.125 * (ny[i  ][j  ][k  ]+ny[i+1][j  ][k  ]
                          +ny[i  ][j+1][k  ]+ny[i+1][j+1][k  ]
                          +ny[i  ][j  ][k+1]+ny[i+1][j  ][k+1]
                          +ny[i  ][j+1][k+1]+ny[i+1][j+1][k+1]);
    real nzcent = 0.125 * (nz[i  ][j  ][k  ]+nz[i+1][j  ][k  ]
                          +nz[i  ][j+1][k  ]+nz[i+1][j+1][k  ]
                          +nz[i  ][j  ][k+1]+nz[i+1][j  ][k+1]
                          +nz[i  ][j+1][k+1]+nz[i+1][j+1][k+1]);
    real nmagcent = sqrt(nxcent*nxcent + nycent*nycent + nzcent*nzcent)+epsnorm;

    /* -div(n) at cell center */
    kappa[i][j][k]= -0.25 * ((nx[i+1][j  ][k  ]-nx[i][j  ][k  ])/dx
                            +(nx[i+1][j+1][k  ]-nx[i][j+1][k  ])/dx
                            +(nx[i+1][j  ][k+1]-nx[i][j  ][k+1])/dx
                            +(nx[i+1][j+1][k+1]-nx[i][j+1][k+1])/dx)
                  - 0.25 * ( (ny[i  ][j+1][k  ]-ny[i  ][j][k  ])/dy
                            +(ny[i+1][j+1][k  ]-ny[i+1][j][k  ])/dy
                            +(ny[i  ][j+1][k+1]-ny[i  ][j][k+1])/dy
                            +(ny[i+1][j+1][k+1]-ny[i+1][j][k+1])/dy)
                  - 0.25 * ( (nz[i  ][j  ][k+1]-nz[i  ][j  ][k])/dz
                            +(nz[i  ][j+1][k+1]-nz[i  ][j+1][k])/dz
                            +(nz[i+1][j  ][k+1]-nz[i+1][j  ][k])/dz
                            +(nz[i+1][j+1][k+1]-nz[i+1][j+1][k])/dz);

    /* (n/|n|.nabla)|n| */
    kappa[i][j][k] += nxcent/nmagcent
                   *0.25 * ( (nmag[i+1][j  ][k  ]-nmag[i][j  ][k  ])/dx
                            +(nmag[i+1][j+1][k  ]-nmag[i][j+1][k  ])/dx
                            +(nmag[i+1][j  ][k+1]-nmag[i][j  ][k+1])/dx
                            +(nmag[i+1][j+1][k+1]-nmag[i][j+1][k+1])/dx)
                   + nycent/nmagcent
                   *0.25 * ( (nmag[i  ][j+1][k  ]-nmag[i  ][j][k  ])/dy
                            +(nmag[i+1][j+1][k  ]-nmag[i+1][j][k  ])/dy
                            +(nmag[i  ][j+1][k+1]-nmag[i  ][j][k+1])/dy
                            +(nmag[i+1][j+1][k+1]-nmag[i+1][j][k+1])/dy)
                   + nzcent/nmagcent
                   *0.25 * ( (nmag[i  ][j  ][k+1]-nmag[i  ][j  ][k])/dz
                            +(nmag[i  ][j+1][k+1]-nmag[i  ][j+1][k])/dz
                            +(nmag[i+1][j  ][k+1]-nmag[i+1][j  ][k])/dz
                            +(nmag[i+1][j+1][k+1]-nmag[i+1][j+1][k])/dz);
    kappa[i][j][k] = kappa[i][j][k]/nmagcent;
  }
#else
  /* cell-center base */
  gradphic();
  for_ijk(i,j,k) {
    real nw = nx[i-1][j][k];
    real ne = nx[i+1][j][k];
    real ns = ny[i][j-1][k];
    real nn = ny[i][j+1][k];
    real nb = nz[i][j][k-1];
    real nt = nz[i][j][k+1];

    kappa[i][j][k]=-((ne-nw)/(phi.dxw(i)+phi.dxe(i))
                    +(nn-ns)/(phi.dys(j)+phi.dyn(j))
                    +(nt-nb)/(phi.dzb(k)+phi.dzt(k)));
  }
#endif

  insert_bc_kappa();
  insert_bc_dist(kappa);
  kappa.exchange();

  return;
}
