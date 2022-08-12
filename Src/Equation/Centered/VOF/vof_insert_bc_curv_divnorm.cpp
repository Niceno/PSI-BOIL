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
void VOF::insert_bc_curv_divnorm() {
/***************************************************************************//**
*  \brief Calculate curvature on wall boundary condition.
*******************************************************************************/

#if 0
  boil::oout<<"VOF::insert_bc_curv_divnorm: obsolete code. Use height-functions instead. "
            <<"Exiting."
            <<boil::endl;
  exit(0);
#endif

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
  // tempflag=1 at "3 <= height-function <=4", which is set in curv_HF.
  if(phi.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=si();
    for_vjk(kappa,j,k) {
      if (tempflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i-1][j][k]=kappa[i][j][k];
      } else {
        tempflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::imax(), BndType::wall())) {
    int i=ei();
    for_vjk(kappa,j,k) {
      if (tempflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i+1][j][k]=kappa[i][j][k];
      } else {
        tempflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::jmin(), BndType::wall())) {
    int j=sj();
    for_vik(kappa,i,k) {
      if (tempflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j-1][k]=kappa[i][j][k];
      } else {
        tempflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::jmax(), BndType::wall())) {
    int j=ej();
    for_vik(kappa,i,k) {
      if (tempflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j+1][k]=kappa[i][j][k];
      } else {
        tempflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::kmin(), BndType::wall())) {
    int k=sk();
    for_vij(kappa,i,j) {
      if (tempflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j][k-1]=kappa[i][j][k];
      } else {
        tempflag[i][j][k]=0;
      }
    }
  }

  if(phi.bc().type_here(Dir::kmax(), BndType::wall())) {
    int k=ek();
    for_vij(kappa,i,j) {
      if (tempflag[i][j][k]==1) {
        kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
        kappa[i][j][k+1]=kappa[i][j][k];
      } else {
        tempflag[i][j][k]=0;
      }
    }
  }

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    if (tempflag[i][j][k]==1) {
      kappa[i][j][k] = divnorm(nx,ny,nz,i,j,k);
    } else {
        tempflag[i][j][k]=0;
    }
  }

  tempflag.exchange();
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
    tempflag2 = tempflag;
    int i=si();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::imin(), BndType::wall())) {
      for_vjk(kappa,j,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (tempflag[i][j][k]==0) {
          int inb =  min(1,tempflag[i][j-1][k]) + min(1,tempflag[i][j+1][k])
                   + min(1,tempflag[i][j][k-1]) + min(1,tempflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,tempflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,tempflag[i][j+1][k])) * kappa[i][j+1][k]
                             + real(min(1,tempflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,tempflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
    }
  }

  /*--------+
  |  i-max  |
  +--------*/
  if(phi.bc().type(Dir::imax(), BndType::wall())) {
    stmp = kappa;
    tempflag2 = tempflag;
    int i=ei();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::imax(), BndType::wall())) {
      for_vjk(kappa,j,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (tempflag[i][j][k]==0) {
          int inb =  min(1,tempflag[i][j-1][k]) + min(1,tempflag[i][j+1][k])
                   + min(1,tempflag[i][j][k-1]) + min(1,tempflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,tempflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,tempflag[i][j+1][k])) * kappa[i][j+1][k]
                             + real(min(1,tempflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,tempflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
    }
  }

  /*--------+
  |  j-min  |
  +--------*/
  if(phi.bc().type(Dir::jmin(), BndType::wall())) {
    stmp = kappa;
    tempflag2 = tempflag;
    int j=sj();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::jmin(), BndType::wall())) {
      for_vik(kappa,i,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (tempflag[i][j][k]==0) {
          int inb =  min(1,tempflag[i-1][j][k]) + min(1,tempflag[i+1][j][k])
                   + min(1,tempflag[i][j][k-1]) + min(1,tempflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,tempflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,tempflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,tempflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,tempflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
    }
  }

  /*--------+
  |  j-max  |
  +--------*/
  if(phi.bc().type(Dir::jmax(), BndType::wall())) {
    stmp = kappa;
    tempflag2 = tempflag;
    int j=ej();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::jmax(), BndType::wall())) {
      for_vik(kappa,i,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if (tempflag[i][j][k]==0) {
          int inb =  min(1,tempflag[i-1][j][k]) + min(1,tempflag[i+1][j][k])
                   + min(1,tempflag[i][j][k-1]) + min(1,tempflag[i][j][k+1]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,tempflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,tempflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,tempflag[i][j][k-1])) * kappa[i][j][k-1]
                             + real(min(1,tempflag[i][j][k+1])) * kappa[i][j][k+1])
                             /real(inb);
              tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
    }
  }

  /*--------+
  |  k-min  |
  +--------*/
  if(phi.bc().type(Dir::kmin(), BndType::wall())) {
    stmp = kappa;
    tempflag2 = tempflag;
    int k=sk();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::kmin(), BndType::wall())) {
      for_vij(kappa,i,j) {
        if(dom->ibody().off(i,j,k)) continue;
        if (tempflag[i][j][k]==0) {
          int inb =  min(1,tempflag[i-1][j][k]) + min(1,tempflag[i+1][j][k])
                   + min(1,tempflag[i][j-1][k]) + min(1,tempflag[i][j+1][k]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,tempflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,tempflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,tempflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,tempflag[i][j+1][k])) * kappa[i][j+1][k])
                             /real(inb);
              tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
    }
  }

  /*--------+
  |  k-max  |
  +--------*/
  if(phi.bc().type(Dir::kmax(), BndType::wall())) {
    stmp = kappa;
    tempflag2 = tempflag;
    int k=ek();
    for(int iloop=1; iloop<=mloop; iloop++) {
      if(phi.bc().type_here(Dir::kmax(), BndType::wall())) {
      for_vij(kappa,i,j) {
        if(dom->ibody().off(i,j,k)) continue;
        if (tempflag[i][j][k]==0) {
          int inb =  min(1,tempflag[i-1][j][k]) + min(1,tempflag[i+1][j][k])
                   + min(1,tempflag[i][j-1][k]) + min(1,tempflag[i][j+1][k]);
          if (inb >= 1) {
              stmp[i][j][k] = (real(min(1,tempflag[i-1][j][k])) * kappa[i-1][j][k]
                             + real(min(1,tempflag[i+1][j][k])) * kappa[i+1][j][k]
                             + real(min(1,tempflag[i][j-1][k])) * kappa[i][j-1][k]
                             + real(min(1,tempflag[i][j+1][k])) * kappa[i][j+1][k])
                             /real(inb);
              tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolated
          }
        }
      }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
    }
  }

  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
  if(dom->ibody().ncall() >= 1) {
    stmp = kappa;
    tempflag2 = tempflag;
    for(int iloop=1; iloop<=mloop; iloop++) {
      for(int cc=0; cc<dom->ibody().nccells(); cc++){
        int i,j,k;
        dom->ibody().ijk(cc,&i,&j,&k);
    /*--------------------------------------+
    |  Note:  (i, j, k) is in fluid domain  |
    +--------------------------------------*/
        if (tempflag[i][j][k]==0) {  // to be extrapolated
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

          int inb = ical*( min(1,tempflag[i-1][j][k]) + min(1,tempflag[i+1][j][k]))
                  + jcal*( min(1,tempflag[i][j-1][k]) + min(1,tempflag[i][j+1][k]))
                  + kcal*( min(1,tempflag[i][j][k-1]) + min(1,tempflag[i][j][k+1]));
          if (inb >= 1) {
            stmp[i][j][k]=(ical*(real(min(1,tempflag[i-1][j][k]))*kappa[i-1][j][k]
                                +real(min(1,tempflag[i+1][j][k]))*kappa[i+1][j][k])
                          +jcal*(real(min(1,tempflag[i][j-1][k]))*kappa[i][j-1][k]
                                +real(min(1,tempflag[i][j+1][k]))*kappa[i][j+1][k])
                          +kcal*(real(min(1,tempflag[i][j][k-1]))*kappa[i][j][k-1]
                               +real(min(1,tempflag[i][j][k+1]))*kappa[i][j][k+1]))
                          /real(inb);
            tempflag2[i][j][k] = 2;  // tempflag=2 for extrapolation
          }
        }
      }
      stmp.exchange();
      tempflag2.exchange();
      kappa = stmp;
      tempflag = tempflag2;
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

