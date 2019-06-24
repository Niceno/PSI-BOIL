#include "vof.h"
#include <iomanip>

void nwall(const Scalar & sca, const real & cangle
         , const real & nwlx, const real & nwly, const real & nwlz
         , const int & i, const int & j, const int & k
         , real nout[]);

real divnorm(const Scalar & nx, const Scalar & ny,const Scalar & nz
            ,const int & i, const int & j, const int & k);

/******************************************************************************/
void VOF::bdcurv() {
/***************************************************************************//**
*  \brief Calculate curvature on wall boundary condition.
*******************************************************************************/

  /* calculate normal vector at vertex */
  gradphi(phi);

  /* cal normal vector on wall considering contact angle */
  wall_norm(phi);

  /* cal normal vector on immersed boundary considering contact angle */
  if(dom->ibody().nccells() > 0) {
    ib_norm(phi);
  }

  /*---------------------------------------------+
  |  calculate curvature. curvature=-div(norm)   |
  +---------------------------------------------*/
  if(phi.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=si();
    for_vjk(kappa,j,k) {
      if (iflag[i][j][k]!=0) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i-1][j][k]=kappa[i][j][k];
      }
    }
  }

  if(phi.bc().type_here(Dir::imax(), BndType::wall())) {
    int i=ei();
    for_vjk(kappa,j,k) {
      if (iflag[i][j][k]!=0) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i+1][j][k]=kappa[i][j][k];
      }
    }
  }

  if(phi.bc().type_here(Dir::jmin(), BndType::wall())) {
    int j=sj();
    for_vik(kappa,i,k) {
      if (iflag[i][j][k]!=0) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j-1][k]=kappa[i][j][k];
      }
    }
  }

  if(phi.bc().type_here(Dir::jmax(), BndType::wall())) {
    int j=ej();
    for_vik(kappa,i,k) {
      if (iflag[i][j][k]!=0) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j+1][k]=kappa[i][j][k];
      }
    }
  }

  if(phi.bc().type_here(Dir::kmin(), BndType::wall())) {
    int k=sk();
    for_vij(kappa,i,j) {
      if (iflag[i][j][k]!=0) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j][k-1]=kappa[i][j][k];
      }
    }
  }

  if(phi.bc().type_here(Dir::kmax(), BndType::wall())) {
    int k=ek();
    for_vij(kappa,i,j) {
      if (iflag[i][j][k]!=0) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j][k+1]=kappa[i][j][k];
      }
    }
  }

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    if (iflag[i][j][k]!=0) {
      kappa[i][j][k] = divnorm(nx,ny,nz,i,j,k);
    }
  }

  kappa.exchange();

  return;
}

/******************************************************************************/
void nwall(const Scalar & sca, const real & cangle
         , const real & nwlx, const real & nwly, const real & nwlz
         , const int & i, const int & j, const int & k
         , real nout[]){
  real epsnorm=1.0e-12;
  real dx=sca.dxw(i);
  real dy=sca.dys(j);
  real dz=sca.dzb(k);
  real dphidx = 0.25*( (sca[i][j  ][k  ]-sca[i-1][j  ][k  ])/dx
                      +(sca[i][j-1][k  ]-sca[i-1][j-1][k  ])/dx
                      +(sca[i][j  ][k-1]-sca[i-1][j  ][k-1])/dx
                      +(sca[i][j-1][k-1]-sca[i-1][j-1][k-1])/dx);
  real dphidy = 0.25*( (sca[i  ][j][k  ]-sca[i  ][j-1][k  ])/dy
                      +(sca[i-1][j][k  ]-sca[i-1][j-1][k  ])/dy
                      +(sca[i  ][j][k-1]-sca[i  ][j-1][k-1])/dy
                      +(sca[i-1][j][k-1]-sca[i-1][j-1][k-1])/dy);
  real dphidz = 0.25*( (sca[i  ][j  ][k]-sca[i  ][j  ][k-1])/dz
                      +(sca[i  ][j-1][k]-sca[i  ][j-1][k-1])/dz
                      +(sca[i-1][j  ][k]-sca[i-1][j  ][k-1])/dz
                      +(sca[i-1][j-1][k]-sca[i-1][j-1][k-1])/dz);
  /* nwl.grad(sca) */
  real tmp = nwlx*dphidx + nwly*dphidy + nwlz*dphidz;
  /* grad(phi)-(nwl.grad(phi))nwl */
  real tmpx= dphidx - tmp*nwlx;
  real tmpy= dphidy - tmp*nwly;
  real tmpz= dphidz - tmp*nwlz;
  real magn= sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz) + epsnorm;
  /* unit vector tangential to the wall */
  real ntanx = tmpx/magn;
  real ntany = tmpy/magn;
  real ntanz = tmpz/magn;
  /* n_wall = nwl * cos(cangle) + ntan * sin(cangle) */
  nout[0] = nwlx*cos(cangle) + ntanx*sin(cangle);
  nout[1] = nwly*cos(cangle) + ntany*sin(cangle);
  nout[2] = nwlz*cos(cangle) + ntanz*sin(cangle);
}

/******************************************************************************/
real divnorm(const Scalar & nx, const Scalar & ny,const Scalar & nz
            ,const int & i, const int & j, const int & k){

  real nw = 0.25* (nx[i  ][j][k  ]+nx[i  ][j+1][k  ]
                  +nx[i  ][j][k+1]+nx[i  ][j+1][k+1]);
  real ne = 0.25* (nx[i+1][j][k  ]+nx[i+1][j+1][k  ]
                  +nx[i+1][j][k+1]+nx[i+1][j+1][k+1]);
  real ns = 0.25* (ny[i][j  ][k  ]+ny[i+1][j  ][k  ]
                  +ny[i][j  ][k+1]+ny[i+1][j  ][k+1]);
  real nn = 0.25* (ny[i][j+1][k  ]+ny[i+1][j+1][k  ]
                  +ny[i][j+1][k+1]+ny[i+1][j+1][k+1]);
  real nb = 0.25* (nz[i][j  ][k  ]+nz[i+1][j  ][k  ]
                  +nz[i][j+1][k  ]+nz[i+1][j+1][k  ]);
  real nt = 0.25* (nz[i][j  ][k+1]+nz[i+1][j  ][k+1]
                  +nz[i][j+1][k+1]+nz[i+1][j+1][k+1]);
  return(-((ne-nw)/(nx.dxc(i))
          +(nn-ns)/(ny.dyc(j))
          +(nt-nb)/(nz.dzc(k))));
}

