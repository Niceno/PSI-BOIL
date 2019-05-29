#include "vof.h"
#include <iomanip>

void nwall(const Scalar & sca, const real & theta
         , const real & nwlx, const real & nwly, const real & nwlz
         , const int & i, const int & j, const int & k
         , real nout[]);

real divnorm(const Scalar & nx, const Scalar & ny,const Scalar & nz
            ,const int & i, const int & j, const int & k);

/******************************************************************************/
void VOF::bdcurv(const Scalar & sca, const real & theta) {
/***************************************************************************//**
*  \brief Calculate curvature on wall boundary condition.
*******************************************************************************/

  gradphi(sca);
  /*---------------------------------------------+
  |  convert grad(phi) to grad(phi)/|grad(phi)|  |
  |   for all node including wall boundary.      |
  +---------------------------------------------*/
  for(int i=1; i<=nx.ei()+1; i++) {
    for(int j=1; j<=nx.ej()+1; j++) {
      for(int k=1; k<=nx.ek()+1; k++) {
        real magn = sqrt(nx[i][j][k]*nx[i][j][k]
                       + ny[i][j][k]*ny[i][j][k]
                       + nz[i][j][k]*nz[i][j][k]) + epsnorm;
        nx[i][j][k] /= magn;
        ny[i][j][k] /= magn;
        nz[i][j][k] /= magn;
      }
    }
  }

  /*----------------------+
  |  wall adhesion        |
  |  calculate n on wall  |
  +----------------------*/
  real nwlx, nwly, nwlz;
  real nout[3];

  /* i-min plane */
  if(sca.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=1;
    for(int j=1; j<=sca.ej()+1; j++){
      for(int k=1; k<=sca.ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx =-1.0;
        nwly = 0.0;
        nwlz = 0.0;
        nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* i-max plane */
  if(sca.bc().type_here(Dir::imax(), BndType::wall())) {
    int i=ei()+1;
    for(int j=1; j<=sca.ej()+1; j++){
      for(int k=1; k<=sca.ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx = 1.0;
        nwly = 0.0;
        nwlz = 0.0;
        nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* j-min plane */
  if(sca.bc().type_here(Dir::jmin(), BndType::wall())) {
    int j=1;
    for(int i=1; i<=sca.ei()+1; i++){
      for(int k=1; k<=sca.ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly =-1.0;
        nwlz = 0.0;
        nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* j-max plane */
  if(sca.bc().type_here(Dir::jmax(), BndType::wall())) {
    int j=ej()+1;
    for(int i=1; i<=sca.ei()+1; i++){
      for(int k=1; k<=sca.ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly = 1.0;
        nwlz = 0.0;
        nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* k-min plane */
  if(sca.bc().type_here(Dir::kmin(), BndType::wall())) {
    int k=1;
    for(int i=1; i<=sca.ei()+1; i++){
      for(int j=1; j<=sca.ej()+1; j++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly = 0.0;
        nwlz =-1.0;
        nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* k-max plane */
  if(sca.bc().type_here(Dir::kmax(), BndType::wall())) {
    int k=ek()+1;
    for(int i=1; i<=sca.ei()+1; i++){
      for(int j=1; j<=sca.ej()+1; j++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly = 0.0;
        nwlz = 1.0;
        nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* line i-min & j-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=1;
    int j=1;
    for(int k=1; k<=sca.ek()+1; k++){
      nwlx =-0.5*sqrt(2.0);
      nwly =-0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-min & j-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=1;
    int j=ej()+1;
    for(int k=1; k<=sca.ek()+1; k++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-min & k-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=1;
    int k=1;
    for(int j=1; j<=sca.ej()+1; j++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-min & k-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=1;
    int k=ek()+1;
    for(int j=1; j<=sca.ej()+1; j++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-max & j-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=1;
    for(int k=1; k<=sca.ek()+1; k++){
      nwlx = 0.5*sqrt(2.0);
      nwly =-0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-max & j-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    for(int k=1; k<=sca.ek()+1; k++){
      nwlx = 0.5*sqrt(2.0);
      nwly = 0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-max & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int k=1;
    for(int j=1; j<=sca.ej()+1; j++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }
 
  /* line i-max & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int k=ek()+1;
    for(int j=1; j<=sca.ej()+1; j++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }
 
  /* line j-min & k-min */
  if(sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=1;
    int k=1;
    for(int i=1; i<=sca.ei()+1; i++){
      nwlx = 0.0;
      nwly =-0.5*sqrt(2.0);
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line j-min & k-max */
  if(sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=1;
    int k=ek()+1;
    for(int i=1; i<=sca.ei()+1; i++){
      nwlx = 0.0;
      nwly =-0.5*sqrt(2.0);
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line j-max & k-min */
  if(sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=ej()+1;
    int k=1;
    for(int i=1; i<=sca.ei()+1; i++){
      nwlx = 0.0;
      nwly = 0.5*sqrt(2.0);
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line j-max & k-max */
  if(sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=ej()+1;
    int k=ek()+1;
    for(int i=1; i<=sca.ei()+1; i++){
      nwlx = 0.0;
      nwly = 0.5*sqrt(2.0);
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* corner i-max & j-min & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=1;
    int k=1;
    nwlx = sqrt(1.0/3.0);
    nwly =-sqrt(1.0/3.0);
    nwlz =-sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-max & j-min & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=1;
    int k=ek()+1;
    nwlx = sqrt(1.0/3.0);
    nwly =-sqrt(1.0/3.0);
    nwlz = sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-min & j-min & k-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=1;
    int j=1;
    int k=1;
    nwlx =-sqrt(1.0/3.0);
    nwly =-sqrt(1.0/3.0);
    nwlz =-sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-min & j-min & k-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=1;
    int j=1;
    int k=ek()+1;
    nwlx =-sqrt(1.0/3.0);
    nwly =-sqrt(1.0/3.0);
    nwlz = sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-min & j-max & k-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=1;
    int j=ej()+1;
    int k=1;
    nwlx =-sqrt(1.0/3.0);
    nwly = sqrt(1.0/3.0);
    nwlz =-sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-min & j-max & k-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=1;
    int j=ej()+1;
    int k=ek()+1;
    nwlx =-sqrt(1.0/3.0);
    nwly = sqrt(1.0/3.0);
    nwlz = sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-max & j-min & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=1;
    int k=1;
    nwlx = sqrt(1.0/3.0);
    nwly =-sqrt(1.0/3.0);
    nwlz =-sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-max & j-min & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=1;
    int k=ek()+1;
    nwlx = sqrt(1.0/3.0);
    nwly =-sqrt(1.0/3.0);
    nwlz = sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-max & j-max & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int k=1;
    nwlx = sqrt(1.0/3.0);
    nwly = sqrt(1.0/3.0);
    nwlz =-sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /* corner i-max & j-max & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int k=ek()+1;
    nwlx = sqrt(1.0/3.0);
    nwly = sqrt(1.0/3.0);
    nwlz = sqrt(1.0/3.0);
    nwall(sca,theta,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  /*---------------------------------------------+
  |  calculate curvature. curvature=-div(norm)   |
  +---------------------------------------------*/
  if(sca.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=1;
    for_vjk(kappa,j,k) {
      kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
      kappa[i-1][j][k]=kappa[i][j][k];
    }
  }

  if(sca.bc().type_here(Dir::imax(), BndType::wall())) {
    int i=ei();
    for_vjk(kappa,j,k) {
      kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
      kappa[i+1][j][k]=kappa[i][j][k];
    }
  }

  if(sca.bc().type_here(Dir::jmin(), BndType::wall())) {
    int j=1;
    for_vik(kappa,i,k) {
      kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
      kappa[i][j-1][k]=kappa[i][j][k];
    }
  }

  if(sca.bc().type_here(Dir::jmax(), BndType::wall())) {
    int j=ej();
    for_vik(kappa,i,k) {
      kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
      kappa[i][j+1][k]=kappa[i][j][k];
    }
  }

  if(sca.bc().type_here(Dir::kmin(), BndType::wall())) {
    int k=1;
    for_vij(kappa,i,j) {
      kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
      kappa[i][j][k-1]=kappa[i][j][k];
    }
  }

  if(sca.bc().type_here(Dir::kmax(), BndType::wall())) {
    int k=ek();
    for_vij(kappa,i,j) {
      kappa[i][j][k]=divnorm(nx,ny,nz,i,j,k);
      kappa[i][j][k+1]=kappa[i][j][k];
    }
  }

  return;
}

/******************************************************************************/
void nwall(const Scalar & sca, const real & theta
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
  /* nwl.grad(phi) */
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
  /* n_wall = nwl * cos(theta) + ntan * sin(theta) */
  nout[0] = nwlx*cos(theta) + ntanx*sin(theta);
  nout[1] = nwly*cos(theta) + ntany*sin(theta);
  nout[2] = nwlz*cos(theta) + ntanz*sin(theta);
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

