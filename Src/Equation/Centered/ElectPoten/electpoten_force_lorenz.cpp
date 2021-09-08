#include "electpoten.h"

//real outerX(int i, int j, int k, const Vector *v1, const Vector *v2);
//real outerY(int i, int j, int k, const Vector *v1, const Vector *v2);
//real outerZ(int i, int j, int k, const Vector *v1, const Vector *v2);
/******************************************************************************/
void ElectPoten::force_lorenz(Vector * xyz) {

  for_m(m){
    for_vmijk((*xyz),m,i,j,k){
      real vect_A[3],vect_B[3],vect_C[3];
      valueFace(&J,m,i,j,k,&vect_A[0],&vect_A[1],&vect_A[2]);
      valueFace(&B,m,i,j,k,&vect_B[0],&vect_B[1],&vect_B[2]);
      boil::crossProduct(vect_C, vect_A, vect_B);
      (*xyz)[m][i][j][k] += vect_C[~m] * (*xyz).dV(m,i,j,k);
    }
  }

  xyz->exchange();

  return;
}

#if 0
/******************************************************************************/
real outerX(int i, int j, int k, const Vector *A, const Vector *B){
  Comp my = Comp::v();
  Comp mz = Comp::w();
  real Ay,Az,By,Bz;
  Ay=0.25*((*A)[my][i-1][j][k]+(*A)[my][i-1][j+1][k]
          +(*A)[my][i  ][j][k]+(*A)[my][i  ][j+1][k]);
  Az=0.25*((*A)[mz][i-1][j][k]+(*A)[mz][i-1][j][k+1]
          +(*A)[mz][i  ][j][k]+(*A)[mz][i  ][j][k+1]);
  By=0.25*((*B)[my][i-1][j][k]+(*B)[my][i-1][j+1][k]
          +(*B)[my][i  ][j][k]+(*B)[my][i  ][j+1][k]);
  Bz=0.25*((*B)[mz][i-1][j][k]+(*B)[mz][i-1][j][k+1]
          +(*B)[mz][i  ][j][k]+(*B)[mz][i  ][j][k+1]);
  return Ay*Bz-Az*By;
}
/******************************************************************************/
real outerY(int i, int j, int k, const Vector *A, const Vector *B){
  Comp mx = Comp::x();
  Comp mz = Comp::w();
  real Ax,Az,Bx,Bz;
  Ax=0.25*((*A)[mx][i][j-1][k]+(*A)[mx][i+1][j-1][k]
          +(*A)[mx][i][j  ][k]+(*A)[mx][i+1][j  ][k]);
  Az=0.25*((*A)[mz][i][j-1][k]+(*A)[mz][i][j-1][k+1]
          +(*A)[mz][i][j  ][k]+(*A)[mz][i][j  ][k+1]);
  Bx=0.25*((*B)[mx][i][j-1][k]+(*B)[mx][i+1][j-1][k]
          +(*B)[mx][i][j  ][k]+(*B)[mx][i+1][j  ][k]);
  Bz=0.25*((*B)[mz][i][j-1][k]+(*B)[mz][i][j-1][k+1]
          +(*B)[mz][i][j  ][k]+(*B)[mz][i][j  ][k+1]);
  return Az*Bx-Ax*Bz;
}
/******************************************************************************/
real outerX(int i, int j, int k, const Vector *A, const Vector *B){
  Comp my = Comp::v();
  Comp mz = Comp::w();
  real Ay,Az,By,Bz;
  Ay=0.25*((*A)[my][i-1][j][k]+(*A)[my][i-1][j+1][k]
          +(*A)[my][i  ][j][k]+(*A)[my][i  ][j+1][k]);
  Az=0.25*((*A)[mz][i-1][j][k]+(*A)[mz][i-1][j][k+1]
          +(*A)[mz][i  ][j][k]+(*A)[mz][i  ][j][k+1]);
  By=0.25*((*B)[my][i-1][j][k]+(*B)[my][i-1][j+1][k]
          +(*B)[my][i  ][j][k]+(*B)[my][i  ][j+1][k]);
  Bz=0.25*((*B)[mz][i-1][j][k]+(*B)[mz][i-1][j][k+1]
          +(*B)[mz][i  ][j][k]+(*B)[mz][i  ][j][k+1]);
  return Ay*Bz-Az*By;
}
#endif

