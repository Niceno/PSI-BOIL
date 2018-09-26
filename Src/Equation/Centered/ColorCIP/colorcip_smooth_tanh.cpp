#include "colorcip.h"
#include <iomanip>

/******************************************************************************/
void ColorCIP::smooth_tanh(const Scalar & sca
                        , Scalar & scb
                        , const int itnum) {
/*----------------------------------+
| smooth sca by hyperbolic tangent  |
|   input:sca, itnum                |
|   output:scb                      |
|   temporary:stmp                  |
+----------------------------------*/

  int iallprd=0;
  /* calculate distance function for boundary condition of tanh-profile */
  if(  sca.bc().type(Dir::imin(), BndType::periodic())
    && sca.bc().type(Dir::imax(), BndType::periodic())
    && sca.bc().type(Dir::jmin(), BndType::periodic())
    && sca.bc().type(Dir::jmax(), BndType::periodic())
    && sca.bc().type(Dir::kmin(), BndType::periodic())
    && sca.bc().type(Dir::kmax(), BndType::periodic())){
    iallprd=1;
  } else {
    smooth_ls(sca, kappa, 8);
  }

  /* set constants */
  const real phimin=-1.0;
  const real phimax= 1.0;
  const real dtau=std::min(4.0*dxmin*dxmin/1.0,ww*ww);

  /* cut off */
  for_aijk(i,j,k)
    scb[i][j][k]=std::max(0.0,(std::min(1.0,sca[i][j][k])));

  /* convert color-function to phi */
  for_aijk(i,j,k)
    scb[i][j][k]=(phimax-phimin)*scb[i][j][k]+phimin;

  /* iterative calculation */
  for(int it=0; it<itnum; it++){

    /* Laplace(phi) */
#if 0
    /* explicit */
    for_ijk(i,j,k){
      stmp[i][j][k] = (scb[i-1][j][k]+scb[i+1][j][k]-2.0*scb[i][j][k])
                      /(sca.dxc(i)*sca.dxc(i))
                     +(scb[i][j-1][k]+scb[i][j+1][k]-2.0*scb[i][j][k])
                      /(sca.dyc(j)*sca.dyc(j))
                     +(scb[i][j][k-1]+scb[i][j][k+1]-2.0*scb[i][j][k])
                      /(sca.dzc(k)*sca.dzc(k));
    }
#else
    /* implicit */
    for_ijk(i,j,k){
      stmp[i][j][k] = (scb[i-1][j][k]+scb[i+1][j][k])
                      /(sca.dxc(i)*sca.dxc(i))
                     +(scb[i][j-1][k]+scb[i][j+1][k])
                      /(sca.dyc(j)*sca.dyc(j))
                     +(scb[i][j][k-1]+scb[i][j][k+1])
                      /(sca.dzc(k)*sca.dzc(k));
    }
#endif
    /* grad(phi) */
    gradphi(scb);

    /* phi*(1-phi**2)/W**2  explicit */
#if 0
    for_ijk(i,j,k)
      stmp[i][j][k] += scb[i][j][k]*(1.0-scb[i][j][k]*scb[i][j][k])/(ww*ww);
#endif

    /* node point */
    for_ijk(i,j,k) {
      real dx=clr.dxc(i);
      real dy=clr.dyc(j);
      real dz=clr.dzc(k);
      /* |grad(phi)| at cell center */
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
      real nmagcent = sqrt(nxcent*nxcent + nycent*nycent + nzcent*nzcent)
                    + epsnorm;

    /* -|grad(phi)|*div(grad(phi)/|grad(phi)|) at cell center */
      stmp[i][j][k] -= nmagcent *
                     (0.25 * ( (nx[i+1][j  ][k  ]/nmag[i+1][j  ][k  ]
                               -nx[i  ][j  ][k  ]/nmag[i  ][j  ][k  ])/dx
                              +(nx[i+1][j+1][k  ]/nmag[i+1][j+1][k  ]
                               -nx[i  ][j+1][k  ]/nmag[i  ][j+1][k  ])/dx
                              +(nx[i+1][j  ][k+1]/nmag[i+1][j  ][k+1]
                               -nx[i  ][j  ][k+1]/nmag[i  ][j  ][k+1])/dx
                              +(nx[i+1][j+1][k+1]/nmag[i+1][j+1][k+1]
                               -nx[i  ][j+1][k+1]/nmag[i  ][j+1][k+1])/dx)
                     +0.25 * ( (ny[i  ][j+1][k  ]/nmag[i  ][j+1][k  ]
                               -ny[i  ][j  ][k  ]/nmag[i  ][j  ][k  ])/dy
                              +(ny[i+1][j+1][k  ]/nmag[i+1][j+1][k  ]
                               -ny[i+1][j  ][k  ]/nmag[i+1][j  ][k  ])/dy
                              +(ny[i  ][j+1][k+1]/nmag[i  ][j+1][k+1]
                               -ny[i  ][j  ][k+1]/nmag[i  ][j  ][k+1])/dy
                              +(ny[i+1][j+1][k+1]/nmag[i+1][j+1][k+1]
                               -ny[i+1][j  ][k+1]/nmag[i+1][j  ][k+1])/dy)
                     +0.25 * ( (nz[i  ][j  ][k+1]/nmag[i  ][j  ][k+1]
                               -nz[i  ][j  ][k  ]/nmag[i  ][j  ][k  ])/dz
                              +(nz[i  ][j+1][k+1]/nmag[i  ][j+1][k+1]
                               -nz[i  ][j+1][k  ]/nmag[i  ][j+1][k  ])/dz
                              +(nz[i+1][j  ][k+1]/nmag[i+1][j  ][k+1]
                               -nz[i+1][j  ][k  ]/nmag[i+1][j  ][k  ])/dz
                              +(nz[i+1][j+1][k+1]/nmag[i+1][j+1][k+1]
                               -nz[i+1][j+1][k  ]/nmag[i+1][j+1][k  ])/dz));
    }

    /* right hand side */
    for_ijk(i,j,k)
      stmp[i][j][k]=scb[i][j][k]+dtau*stmp[i][j][k];

    /* update */
    for_ijk(i,j,k){
      real diag=1.0
               -dtau* ( (1.0-scb[i][j][k]*scb[i][j][k])/(ww*ww)
                       +(-2.0/(sca.dxc(i)*sca.dxc(i))
                         -2.0/(sca.dyc(j)*sca.dyc(j))
                         -2.0/(sca.dzc(k)*sca.dzc(k))) );
      scb[i][j][k] = stmp[i][j][k]/diag;
    }

    /* set boundary value which is calculated by levelset */
    if(iallprd!=1) insert_bc_tanh(scb,kappa);
    /* set boundary condition */
    insert_bc_ls(scb);
    scb.exchange_all();

  }

  /* convert phi to color-function */
  for_ijk(i,j,k)
    scb[i][j][k]=(scb[i][j][k]-phimin)/(phimax-phimin);

  /* cut off */
  for_ijk(i,j,k)
    scb[i][j][k]=std::max(0.0,(std::min(1.0,scb[i][j][k])));

  /* set boundary condition */
  insert_bc(scb);
  scb.exchange_all();

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_smooth_tanh.cpp,v 1.2 2009/11/12 12:15:48 sato Exp $'/
+-----------------------------------------------------------------------------*/
