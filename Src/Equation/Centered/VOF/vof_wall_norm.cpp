#include "vof.h"
#include <iomanip>

/******************************************************************************/
void VOF::wall_norm(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector on wall
*  Eq. 53 in Brackbill's paper, JCP 1991
*******************************************************************************/

  /*----------------------+
  |  wall adhesion        |
  |  calculate n on wall  |
  +----------------------*/
  real nwlx, nwly, nwlz;
  real nout[3];

  /* i-min plane */
  if(sca.bc().type_here(Dir::imin(), BndType::wall())) {
    int i=si();
    for(int j=sj(); j<=ej()+1; j++){
      for(int k=sk(); k<=ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx =-1.0;
        nwly = 0.0;
        nwlz = 0.0;
        nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* i-max plane */
  if(sca.bc().type_here(Dir::imax(), BndType::wall())) {
    int i=ei()+1;
    for(int j=sj(); j<=ej()+1; j++){
      for(int k=sk(); k<=ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx = 1.0;
        nwly = 0.0;
        nwlz = 0.0;
        nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* j-min plane */
  if(sca.bc().type_here(Dir::jmin(), BndType::wall())) {
    int j=sj();
    for(int i=si(); i<=ei()+1; i++){
      for(int k=sk(); k<=ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly =-1.0;
        nwlz = 0.0;
        nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* j-max plane */
  if(sca.bc().type_here(Dir::jmax(), BndType::wall())) {
    int j=ej()+1;
    for(int i=si(); i<=ei()+1; i++){
      for(int k=sk(); k<=ek()+1; k++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly = 1.0;
        nwlz = 0.0;
        nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* k-min plane */
  if(sca.bc().type_here(Dir::kmin(), BndType::wall())) {
    int k=sk();
    for(int i=si(); i<=ei()+1; i++){
      for(int j=sj(); j<=ej()+1; j++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly = 0.0;
        nwlz =-1.0;
        nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
        nx[i][j][k]=nout[0];
        ny[i][j][k]=nout[1];
        nz[i][j][k]=nout[2];
      }
    }
  }

  /* k-max plane */
  if(sca.bc().type_here(Dir::kmax(), BndType::wall())) {
    int k=ek()+1;
    for(int i=si(); i<=ei()+1; i++){
      for(int j=sj(); j<=ej()+1; j++){
        /* unit normal vector directed into the wall */
        nwlx = 0.0;
        nwly = 0.0;
        nwlz = 1.0;
        nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int k=sk(); k<=ek()+1; k++){
      nwlx =-0.5*sqrt(2.0);
      nwly =-0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int k=sk(); k<=ek()+1; k++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int j=sj(); j<=ej()+1; j++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int j=sj(); j<=ej()+1; j++){
      nwlx =-0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int k=sk(); k<=ek()+1; k++){
      nwlx = 0.5*sqrt(2.0);
      nwly =-0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int k=sk(); k<=ek()+1; k++){
      nwlx = 0.5*sqrt(2.0);
      nwly = 0.5*sqrt(2.0);
      nwlz = 0.0;
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
      nx[i][j][k]=nout[0];
      ny[i][j][k]=nout[1];
      nz[i][j][k]=nout[2];
    }
  }

  /* line i-max & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int k=sk();
    for(int j=sj(); j<=ej()+1; j++){
      nwlx = 0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int j=sj(); j<=ej()+1; j++){
      nwlx = 0.5*sqrt(2.0);
      nwly = 0.0;
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int i=si(); i<=ei()+1; i++){
      nwlx = 0.0;
      nwly =-0.5*sqrt(2.0);
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int i=si(); i<=ei()+1; i++){
      nwlx = 0.0;
      nwly =-0.5*sqrt(2.0);
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int i=si(); i<=ei()+1; i++){
      nwlx = 0.0;
      nwly = 0.5*sqrt(2.0);
      nwlz =-0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    for(int i=si(); i<=ei()+1; i++){
      nwlx = 0.0;
      nwly = 0.5*sqrt(2.0);
      nwlz = 0.5*sqrt(2.0);
      nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
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
    nwall(sca,nwlx,nwly,nwlz,i,j,k,nout);
    nx[i][j][k]=nout[0];
    ny[i][j][k]=nout[1];
    nz[i][j][k]=nout[2];
  }

  return;
}

/******************************************************************************/
void VOF::nwall(const Scalar & sca
         , const real & nwlx, const real & nwly, const real & nwlz
         , const int & i, const int & j, const int & k
         , real nout[]){

  /* dphidx,dphidy,dphidz */
  int ii=i;
  if(i==si() 
     && sca.bc().type_here(Dir::imin(), BndType::wall())){
    ii=si()+1;
  } else if(i==ei()+1
     && sca.bc().type_here(Dir::imax(), BndType::wall())){
    ii=ei();
  }
  int jj=j;
  if(j==sj()
     && sca.bc().type_here(Dir::jmin(), BndType::wall())){
    jj=sj()+1;
  } else if(j==ej()+1
     && sca.bc().type_here(Dir::jmax(), BndType::wall())){
    jj=ej();
  }
  int kk=k;
  if(k==sk()
     && sca.bc().type_here(Dir::kmin(), BndType::wall())){
    kk=sk()+1;
  } else if(k==ek()+1
     && sca.bc().type_here(Dir::kmax(), BndType::wall())){
    kk=ek();
  }

  real dx=sca.dxw(ii);
  real dy=sca.dys(jj);
  real dz=sca.dzb(kk);

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

  if(i==si()
     && sca.bc().type_here(Dir::imin(), BndType::symmetry())){
    dphidx=0.0;
  } else if(i==ei()+1
     && sca.bc().type_here(Dir::imax(), BndType::symmetry())){
    dphidx=0.0;
  }
  if(j==sj()
     && sca.bc().type_here(Dir::jmin(), BndType::symmetry())){
    dphidy=0.0;
  } else if(j==ej()+1
     && sca.bc().type_here(Dir::jmax(), BndType::symmetry())){
    dphidy=0.0;
  }
  if(k==sk()
     && sca.bc().type_here(Dir::kmin(), BndType::symmetry())){
    dphidz=0.0;
  } else if(k==ek()+1
     && sca.bc().type_here(Dir::kmax(), BndType::symmetry())){
    dphidz=0.0;
  }

  /* nwl.grad(phi) */
  real pin = nwlx*dphidx + nwly*dphidy + nwlz*dphidz;
  /* grad(phi)-(nwl.grad(phi))nwl */
  real ntanx= dphidx - pin*nwlx;
  real ntany= dphidy - pin*nwly;
  real ntanz= dphidz - pin*nwlz;
  normalize(ntanx,ntany,ntanz);
  /* n_wall = nwl * cos(cangle) + ntan * sin(cangle) */
  real nnx = nwlx*cos(cangle) + ntanx*sin(cangle);
  real nny = nwly*cos(cangle) + ntany*sin(cangle);
  real nnz = nwlz*cos(cangle) + ntanz*sin(cangle);
  normalize(nnx,nny,nnz);
  nout[0] = nnx;
  nout[1] = nny;
  nout[2] = nnz;
}
