#include "marching_cube.h"
#include <iostream>

/* interface to call marching_cube_surfval using cell coordinates,
   direction and orientation, sig is neg() - negative direction and vice versa */
real MarchingCube::surface(const Sign sig, const Comp & mcomp,
                           const int i, const int j, const int k) {
  CELLFACE face;
  int isum(0);
  real surf;

  if (mcomp == Comp::i()) {
    int ii = i;
    surf = (*clr).dSx(sig,i,j,k);
    if (sig == Sign::neg()) { ii -= 1; }
    for (int m = 0; m != 4; m++) {
      int jj, kk;
      if(m==0)      {jj=j-1; kk=k  ;}
      else if(m==1) {jj=j-1; kk=k-1;}
      else if(m==2) {jj=j  ; kk=k-1;}
      else if(m==3) {jj=j  ; kk=k  ;}
 
      face.val[m] = ((*clr)[ii][jj  ][kk  ]+(*clr)[ii+1][jj  ][kk  ]
                    +(*clr)[ii][jj+1][kk  ]+(*clr)[ii+1][jj+1][kk  ]
                    +(*clr)[ii][jj  ][kk+1]+(*clr)[ii+1][jj  ][kk+1]
                    +(*clr)[ii][jj+1][kk+1]+(*clr)[ii+1][jj+1][kk+1])/8.0;

      if (fabs(face.val[m]-clrsurf) < 1.0e-11) {
        face.val[m]=clrsurf+1.0e-11;
      }
     
      if (face.val[m] > clrsurf) isum++;

      /* x: y->x, z->y */
      face.p[m].x = (*clr).yn(jj+1);
      face.p[m].y = (*clr).zn(kk+1);

    }
  }

  else if (mcomp == Comp::j()) {
    int jj = j;
    surf = (*clr).dSy(sig,i,j,k);
    if (sig == Sign::neg()) { jj -= 1; }
    for (int m = 0; m != 4; m++) {
      int ii, kk;
      if(m==0)      {ii=i-1; kk=k  ;}       
      else if(m==1) {ii=i-1; kk=k-1;}
      else if(m==2) {ii=i  ; kk=k-1;}
      else if(m==3) {ii=i  ; kk=k  ;}

      face.val[m] = ((*clr)[ii][jj  ][kk  ]+(*clr)[ii+1][jj  ][kk  ]
                    +(*clr)[ii][jj+1][kk  ]+(*clr)[ii+1][jj+1][kk  ]
                    +(*clr)[ii][jj  ][kk+1]+(*clr)[ii+1][jj  ][kk+1]
                    +(*clr)[ii][jj+1][kk+1]+(*clr)[ii+1][jj+1][kk+1])/8.0;

      if (fabs(face.val[m]-clrsurf) < 1.0e-11) {
        face.val[m]=clrsurf+1.0e-11; 
      }

      if (face.val[m] > clrsurf) isum++;

      /* y: x->x, z->y */
      face.p[m].x = (*clr).xn(ii+1);
      face.p[m].y = (*clr).zn(kk+1);
    }
  }

  else if (mcomp == Comp::k()) {
    int kk = k;
    surf = (*clr).dSz(sig,i,j,k);
    if (sig == Sign::neg()) { kk -= 1; }
    for (int m = 0; m != 4; m++) {
      int ii, jj;
      if(m==0)      {ii=i-1; jj=j  ;}      
      else if(m==1) {ii=i-1; jj=j-1;}
      else if(m==2) {ii=i  ; jj=j-1;}
      else if(m==3) {ii=i  ; jj=j  ;}

      face.val[m] = ((*clr)[ii][jj  ][kk  ]+(*clr)[ii+1][jj  ][kk  ]
                    +(*clr)[ii][jj+1][kk  ]+(*clr)[ii+1][jj+1][kk  ]
                    +(*clr)[ii][jj  ][kk+1]+(*clr)[ii+1][jj  ][kk+1]
                    +(*clr)[ii][jj+1][kk+1]+(*clr)[ii+1][jj+1][kk+1])/8.0;

      if (fabs(face.val[m]-clrsurf) < 1.0e-11) {
        face.val[m]=clrsurf+1.0e-11; 
      }

      if (face.val[m] > clrsurf) isum++;

      /* z: x->x, y->y */
      face.p[m].x = (*clr).xn(ii+1);
      face.p[m].y = (*clr).yn(jj+1);
    }
  }
 
  real area(0.0);

  if (isum == 4) area += 1.0;
  else if (isum > 0) 
    area += surfval(face, clrsurf)/surf;

  return (area);
}

