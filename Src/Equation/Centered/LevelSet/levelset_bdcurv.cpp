#include "levelset.h"
#include <iomanip>

void nwall(const Scalar & sca, const real & theta
         , const real & nwlx, const real & nwly, const real & nwlz
         , const int & i, const int & j, const int & k
         , real nout[]);

real divnorm(const Scalar & nx, const Scalar & ny,const Scalar & nz
            ,const int & i, const int & j, const int & k);

/******************************************************************************/
void LevelSet::bdcurv(const Scalar & sca, const real & theta) {
/***************************************************************************//**
*  \brief Calculate curvature on wall boundary condition.
*******************************************************************************/

  gradphi();

  /*----------------------+
  |  wall adhesion        |
  |  calculate n on wall  |
  +----------------------*/
  real nwlx, nwly, nwlz;
  real nout[3];

  /* i-min plane */
  if(sca.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=si();
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
    int j=sj();
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
    int k=sk();
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
    int i=si();
    int j=sj();
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
    int i=si();
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
    int i=si();
    int k=sk();
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
    int i=si();
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
    int j=sj();
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
    int k=si();
    for(int j=1; j<=sca.ej()+1; j++){
      nwlx = 0.5*sqrt(2.0);
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
      nwlx = 0.5*sqrt(2.0);
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
    int j=sj();
    int k=sk();
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
    int j=sj();
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
    int k=sk();
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
    int j=sj();
    int k=sk();
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
    int j=sj();
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
    int i=si();
    int j=sj();
    int k=sk();
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
    int i=si();
    int j=sj();
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
    int i=si();
    int j=ej()+1;
    int k=sk();
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
    int i=si();
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
    int j=sj();
    int k=sk();
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
    int j=sj();
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
    int k=sk();
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
    int ii=si();
    for_vjk(kappa,j,k) {
      kappa[ii][j][k]=divnorm(nx,ny,nz,ii,j,k);
      kappa[ii-1][j][k]=kappa[ii][j][k];
    }
  }

  if(sca.bc().type_here(Dir::imax(), BndType::wall())) {
    int ii=ei();
    for_vjk(kappa,j,k) {
      kappa[ii][j][k]=divnorm(nx,ny,nz,ii,j,k);
      kappa[ii+1][j][k]=kappa[ii][j][k];
    }
  }

  if(sca.bc().type_here(Dir::jmin(), BndType::wall())) {
    int jj=sj();
    for_vik(kappa,i,k) {
      kappa[i][jj][k]=divnorm(nx,ny,nz,i,jj,k);
      kappa[i][jj-1][k]=kappa[i][jj][k];
    }
  }

  if(sca.bc().type_here(Dir::jmax(), BndType::wall())) {
    int jj=ej();
    for_vik(kappa,i,k) {
      kappa[i][jj][k]=divnorm(nx,ny,nz,i,jj,k);
      kappa[i][jj+1][k]=kappa[i][jj][k];
    }
  }

  if(sca.bc().type_here(Dir::kmin(), BndType::wall())) {
    int kk=sk();
    for_vij(kappa,i,j) {
      kappa[i][j][kk]=divnorm(nx,ny,nz,i,j,kk);
      kappa[i][j][kk-1]=kappa[i][j][kk];
    }
  }

  if(sca.bc().type_here(Dir::kmax(), BndType::wall())) {
    int kk=ek();
    for_vij(kappa,i,j) {
      kappa[i][j][kk]=divnorm(nx,ny,nz,i,j,kk);
      kappa[i][j][kk+1]=kappa[i][j][kk];
    }
  }

  //std::cout<<"bdcurv:kappa="<<kappa[0][16][16]<<" "<<kappa[1][16][16]<<"\n";
  //std::cout<<"bdcurv:kappa="<<kappa[32][16][16]<<" "<<kappa[33][16][16]<<"\n";

  return;
}

/******************************************************************************/
void nwall(const Scalar & sca, const real & theta
         , const real & nwlx, const real & nwly, const real & nwlz
         , const int & i, const int & j, const int & k
         , real nout[]){

  real epsnorm=1.0e-12;

  /* dphidx,dphidy,dphidz */
  int ii=i;
  if(i==sca.si() 
     && sca.bc().type_here(Dir::imin(), BndType::wall())){
    ii=sca.si()+1;
  } else if(i==sca.ei()+1
     && sca.bc().type_here(Dir::imax(), BndType::wall())){
    ii=sca.ei();
  }
  int jj=j;
  if(j==sca.sj()
     && sca.bc().type_here(Dir::jmin(), BndType::wall())){
    jj=sca.sj()+1;
  } else if(j==sca.ej()+1
     && sca.bc().type_here(Dir::jmax(), BndType::wall())){
    jj=sca.ej();
  }
  int kk=k;
  if(k==sca.sk()
     && sca.bc().type_here(Dir::kmin(), BndType::wall())){
    kk=sca.sk()+1;
  } else if(k==sca.ek()+1
     && sca.bc().type_here(Dir::kmax(), BndType::wall())){
    kk=sca.ek();
  }

  real dx=sca.dxw(ii);
  real dy=sca.dys(jj);
  real dz=sca.dzb(kk);
#if 0
  real dphidx = 0.25*( (sca[ii][jj  ][kk  ]-sca[ii-1][jj  ][kk  ])/dx
                      +(sca[ii][jj-1][kk  ]-sca[ii-1][jj-1][kk  ])/dx
                      +(sca[ii][jj  ][kk-1]-sca[ii-1][jj  ][kk-1])/dx
                      +(sca[ii][jj-1][kk-1]-sca[ii-1][jj-1][kk-1])/dx);
  real dphidy = 0.25*( (sca[ii  ][jj][kk  ]-sca[ii  ][jj-1][kk  ])/dy
                      +(sca[ii-1][jj][kk  ]-sca[ii-1][jj-1][kk  ])/dy
                      +(sca[ii  ][jj][kk-1]-sca[ii  ][jj-1][kk-1])/dy
                      +(sca[ii-1][jj][kk-1]-sca[ii-1][jj-1][kk-1])/dy);
  real dphidz = 0.25*( (sca[ii  ][jj  ][kk]-sca[ii  ][jj  ][kk-1])/dz
                      +(sca[ii  ][jj-1][kk]-sca[ii  ][jj-1][kk-1])/dz
                      +(sca[ii-1][jj  ][kk]-sca[ii-1][jj  ][kk-1])/dz
                      +(sca[ii-1][jj-1][kk]-sca[ii-1][jj-1][kk-1])/dz);
#else
  real dphidx = 0.25*( (sca[ii][j  ][k  ]-sca[ii-1][j  ][k  ])/dx
                      +(sca[ii][j-1][k  ]-sca[ii-1][j-1][k  ])/dx
                      +(sca[ii][j  ][k-1]-sca[ii-1][j  ][k-1])/dx
                      +(sca[ii][j-1][k-1]-sca[ii-1][j-1][k-1])/dx);
  real dphidy = 0.25*( (sca[i  ][jj][k  ]-sca[i  ][jj-1][k  ])/dy
                      +(sca[i-1][jj][k  ]-sca[i-1][jj-1][k  ])/dy
                      +(sca[i  ][jj][k-1]-sca[i  ][jj-1][k-1])/dy
                      +(sca[i-1][jj][k-1]-sca[i-1][jj-1][k-1])/dy);
  real dphidz = 0.25*( (sca[i  ][j  ][kk]-sca[i  ][j  ][kk-1])/dz
                      +(sca[i  ][j-1][kk]-sca[i  ][j-1][kk-1])/dz
                      +(sca[i-1][j  ][kk]-sca[i-1][j  ][kk-1])/dz
                      +(sca[i-1][j-1][kk]-sca[i-1][j-1][kk-1])/dz);
#endif
#if 0
  if(i==2&&j==1&&k==1){
    std::cout<<sca[ii][j  ][k  ]-sca[ii-1][j  ][k  ]<<" "
    <<sca[ii][j-1][k  ]-sca[ii-1][j-1][k  ]<<" "
    <<sca[ii][j  ][k-1]-sca[ii-1][j  ][k-1]<<" "
    <<sca[ii][j-1][k-1]-sca[ii-1][j-1][k-1]<<"\n";
    std::cout<<ii<<" "<<jj<<" "<<kk<<"\n";
  }
#endif

  if(i==sca.si()
     && sca.bc().type_here(Dir::imin(), BndType::symmetry())){
    dphidx=0.0;
  } else if(i==sca.ei()+1
     && sca.bc().type_here(Dir::imax(), BndType::symmetry())){
    dphidx=0.0;
  }
  if(j==sca.sj()
     && sca.bc().type_here(Dir::jmin(), BndType::symmetry())){
    dphidy=0.0;
  } else if(j==sca.ej()+1
     && sca.bc().type_here(Dir::jmax(), BndType::symmetry())){
    dphidy=0.0;
  }
  if(k==sca.sk()
     && sca.bc().type_here(Dir::kmin(), BndType::symmetry())){
    dphidz=0.0;
  } else if(k==sca.ek()+1
     && sca.bc().type_here(Dir::kmax(), BndType::symmetry())){
    dphidz=0.0;
  }

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
#if 0
  if(i==1&&j==17&&k==17){
     std::cout<<"bdcurv "<<ii<<" "<<jj<<" "<<kk<<" "
                   <<dphidx<<" "<<dphidy<<" "<<dphidz<<" "
                   <<nout[0]<<" "<<nout[1]<<" "<<nout[2]<<"\n";
    std::cout<<"nwall "<<nwlx<<" "<<nwly<<" "<<nwlz<<"\n";
  }
#endif
#if 0
  if(i==33&&j==17&&k==17){
     std::cout<<"bdcurv "<<ii<<" "<<jj<<" "<<kk<<" "
                   <<dphidx<<" "<<dphidy<<" "<<dphidz<<" "
                   <<nout[0]<<" "<<nout[1]<<" "<<nout[2]<<"\n";
    std::cout<<"nwall "<<nwlx<<" "<<nwly<<" "<<nwlz<<"\n";
  }
#endif
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

#if 0
  if(i==1&&j==16&&k==16) std::cout<<"bdcurv "
                   <<nw<<" "<<ne<<" "<<ns<<" "<<nn<<" "<<nb<<" "<<nt<<" "
                   <<nx.dxc(i)<<" "<<ny.dyc(j)<<" "<<nz.dzc(k)<<"\n";
#endif
#if 0
  if(i==32&&j==16&&k==16) std::cout<<"bdcurv "
                   <<nw<<" "<<ne<<" "<<ns<<" "<<nn<<" "<<nb<<" "<<nt<<" "
                   <<nx.dxc(i)<<" "<<ny.dyc(j)<<" "<<nz.dzc(k)<<"\n";
#endif
  return(-((ne-nw)/(nx.dxc(i))
          +(nn-ns)/(ny.dyc(j))
          +(nt-nb)/(nz.dzc(k))));

}
