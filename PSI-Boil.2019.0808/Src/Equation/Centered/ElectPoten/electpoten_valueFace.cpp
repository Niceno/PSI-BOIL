#include "electpoten.h"

/******************************************************************************/
void ElectPoten::valueFace(const Vector * A,const Comp m,
               const int i, const int j, const int k, 
               real *ax,real *ay,real *az){
  Comp mi = Comp::u();
  Comp mj = Comp::v();
  Comp mk = Comp::w();
  if(m==Comp::u()){
    *ax = (*A)[mi][i][j][k];
    *ay = 0.25*((*A)[mj][i-1][j][k]+(*A)[mj][i-1][j+1][k]
               +(*A)[mj][i  ][j][k]+(*A)[mj][i  ][j+1][k]);
    *az = 0.25*((*A)[mk][i-1][j][k]+(*A)[mk][i-1][j][k+1]
               +(*A)[mk][i  ][j][k]+(*A)[mk][i  ][j][k+1]);
  } else if(m==Comp::v()){
    *ax = 0.25*((*A)[mi][i][j-1][k]+(*A)[mi][i+1][j-1][k]
               +(*A)[mi][i][j  ][k]+(*A)[mi][i+1][j  ][k]);
    *ay = (*A)[mj][i][j][k];
    *az = 0.25*((*A)[mk][i][j-1][k]+(*A)[mk][i][j-1][k+1]
               +(*A)[mk][i][j  ][k]+(*A)[mk][i][j  ][k+1]);
  } else if(m==Comp::w()){
    *ax = 0.25*((*A)[mi][i][j][k-1]+(*A)[mi][i+1][j][k-1]
               +(*A)[mi][i][j][k  ]+(*A)[mi][i+1][j][k  ]);
    *ay = 0.25*((*A)[mj][i][j][k-1]+(*A)[mj][i][j+1][k-1]
               +(*A)[mj][i][j][k  ]+(*A)[mj][i][j+1][k]);
    *az = (*A)[mk][i][j][k];
  }
}
