#include "vof.h"
#include <iomanip>
using namespace std;

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

  boil::oout<<"VOF::bdcurv: obsolete code. Use height-functions instead. "
            <<"Exiting."
            <<boil::endl;
  exit(0);

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
  // iflag=1 at "3 <= height-function <=4", which is set in curv_HF.
  if(phi.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=si();
    for_vjk(kappa,j,k) {
      if (iflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i-1][j][k]=kappa[i][j][k];
      } else {
        iflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::imax(), BndType::wall())) {
    int i=ei();
    for_vjk(kappa,j,k) {
      if (iflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i+1][j][k]=kappa[i][j][k];
      } else {
        iflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::jmin(), BndType::wall())) {
    int j=sj();
    for_vik(kappa,i,k) {
      if (iflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j-1][k]=kappa[i][j][k];
      } else {
        iflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::jmax(), BndType::wall())) {
    int j=ej();
    for_vik(kappa,i,k) {
      if (iflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j+1][k]=kappa[i][j][k];
      } else {
        iflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::kmin(), BndType::wall())) {
    int k=sk();
    for_vij(kappa,i,j) {
      if (iflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j][k-1]=kappa[i][j][k];
      } else {
        iflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::kmax(), BndType::wall())) {
    int k=ek();
    for_vij(kappa,i,j) {
      if (iflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j][k+1]=kappa[i][j][k];
      } else {
        iflag[i][j][k]=0;
      }
    }
  }

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    if (iflag[i][j][k]==1) {
      kappa[i][j][k] = divnorm(nx,ny,nz,i,j,k);
    } else {
        iflag[i][j][k]=0;
    }
  }

  iflag.exchange();
  kappa.exchange();

  /*---------------------------------+
  |  extrapolate curvature on wall   |
  +---------------------------------*/
  const int mloop = 4;  // number of loop for extrapolation

  /*--------+
  |  i-min  |
  +--------*/
  if(phi.bc().type(Dir::imin(), BndType::wall())) {
    stmp = kappa;
    jflag = iflag;
    int i=si();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::imin(), BndType::wall())) {
      for_vjk(kappa,j,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (iflag[i][j][k]==0) {
          int inb =  min(1,iflag[i][j-1][k]) + min(1,iflag[i][j+1][k])
                   + min(1,iflag[i][j][k-1]) + min(1,iflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,iflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,iflag[i][j+1][k])) * kappa[i][j+1][k]
                             + real(min(1,iflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,iflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              jflag[i][j][k] = 2;  // iflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

  /*--------+
  |  i-max  |
  +--------*/
  if(phi.bc().type(Dir::imax(), BndType::wall())) {
    stmp = kappa;
    jflag = iflag;
    int i=ei();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::imax(), BndType::wall())) {
      for_vjk(kappa,j,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (iflag[i][j][k]==0) {
          int inb =  min(1,iflag[i][j-1][k]) + min(1,iflag[i][j+1][k])
                   + min(1,iflag[i][j][k-1]) + min(1,iflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,iflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,iflag[i][j+1][k])) * kappa[i][j+1][k]
                             + real(min(1,iflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,iflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              jflag[i][j][k] = 2;  // iflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

  /*--------+
  |  j-min  |
  +--------*/
  if(phi.bc().type(Dir::jmin(), BndType::wall())) {
    stmp = kappa;
    jflag = iflag;
    int j=sj();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::jmin(), BndType::wall())) {
      for_vik(kappa,i,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (iflag[i][j][k]==0) {
          int inb =  min(1,iflag[i-1][j][k]) + min(1,iflag[i+1][j][k])
                   + min(1,iflag[i][j][k-1]) + min(1,iflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,iflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,iflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,iflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,iflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              jflag[i][j][k] = 2;  // iflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

  /*--------+
  |  j-max  |
  +--------*/
  if(phi.bc().type(Dir::jmax(), BndType::wall())) {
    stmp = kappa;
    jflag = iflag;
    int j=ej();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::jmax(), BndType::wall())) {
      for_vik(kappa,i,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (iflag[i][j][k]==0) {
          int inb =  min(1,iflag[i-1][j][k]) + min(1,iflag[i+1][j][k])
                   + min(1,iflag[i][j][k-1]) + min(1,iflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,iflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,iflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,iflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,iflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              jflag[i][j][k] = 2;  // iflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

  /*--------+
  |  k-min  |
  +--------*/
  if(phi.bc().type(Dir::kmin(), BndType::wall())) {
    stmp = kappa;
    jflag = iflag;
    int k=sk();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::kmin(), BndType::wall())) {
      for_vij(kappa,i,j) {
        if(dom->ibody().off(i,j,k)) continue;
        if (iflag[i][j][k]==0) {
          int inb =  min(1,iflag[i-1][j][k]) + min(1,iflag[i+1][j][k])
                   + min(1,iflag[i][j-1][k]) + min(1,iflag[i][j+1][k]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,iflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,iflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,iflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,iflag[i][j+1][k])) * kappa[i][j+1][k])
                             /real(inb);
              jflag[i][j][k] = 2;  // iflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

  /*--------+
  |  k-max  |
  +--------*/
  if(phi.bc().type(Dir::kmax(), BndType::wall())) {
    stmp = kappa;
    jflag = iflag;
    int k=ek();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::kmax(), BndType::wall())) {
      for_vij(kappa,i,j) {
        if(dom->ibody().off(i,j,k)) continue;
        if (iflag[i][j][k]==0) {
          int inb =  min(1,iflag[i-1][j][k]) + min(1,iflag[i+1][j][k])
                   + min(1,iflag[i][j-1][k]) + min(1,iflag[i][j+1][k]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,iflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,iflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,iflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,iflag[i][j+1][k])) * kappa[i][j+1][k])
                             /real(inb);
              jflag[i][j][k] = 2;  // iflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
  if(dom->ibody().ncall() >= 1) {
    stmp = kappa;
    jflag = iflag;
    for(int iloop=1; iloop<=mloop; iloop++) {
      for(int cc=0; cc<dom->ibody().nccells(); cc++){
        int i,j,k;
        dom->ibody().ijk(cc,&i,&j,&k);
    /*--------------------------------------+
    |  Note:  (i, j, k) is in fluid domain  |
    +--------------------------------------*/
        if (iflag[i][j][k]==0) {  // to be extrapolated
          /* set direction */    // crude code!!!
          real ux=fabs(dom->ibody().nwx(i,j,k));
          real uy=fabs(dom->ibody().nwy(i,j,k));
          real uz=fabs(dom->ibody().nwz(i,j,k));

          int dirMax=0;
          if(ux>uy && ux>uz) {
            dirMax=0;
          } else if(uy>ux && uy>uz) {
            dirMax=1;
          } else {
            dirMax=2;
          }

          int ical,jcal,kcal;
          ical=jcal=kcal=1;
          if(dirMax==0)      { ical=0; }
          else if(dirMax==1) { jcal=0; }
          else               { kcal=0; }

          int inb = ical*( min(1,iflag[i-1][j][k]) + min(1,iflag[i+1][j][k]))
                  + jcal*( min(1,iflag[i][j-1][k]) + min(1,iflag[i][j+1][k]))
                  + kcal*( min(1,iflag[i][j][k-1]) + min(1,iflag[i][j][k+1]));
          if (inb >= 1) {
            stmp[i][j][k]=(ical*(real(min(1,iflag[i-1][j][k]))*kappa[i-1][j][k]
                                +real(min(1,iflag[i+1][j][k]))*kappa[i+1][j][k])
                          +jcal*(real(min(1,iflag[i][j-1][k]))*kappa[i][j-1][k]
                                +real(min(1,iflag[i][j+1][k]))*kappa[i][j+1][k])
                          +kcal*(real(min(1,iflag[i][j][k-1]))*kappa[i][j][k-1]
                               +real(min(1,iflag[i][j][k+1]))*kappa[i][j][k+1]))
                          /real(inb);
            jflag[i][j][k] = 2;  // iflag=2 for extrapolation
          }
        }
      }
      stmp.exchange();
      jflag.exchange();
      kappa = stmp;
      iflag = jflag;
    }
  }

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

