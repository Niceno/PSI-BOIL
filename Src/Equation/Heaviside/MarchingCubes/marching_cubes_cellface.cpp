#include "marching_cubes.h"
#include <iostream>

/* interface to call heaviside_surface_covered using cell coordinates,
   direction and orientation, sig is neg() - negative direction and vice versa,
   returning area density
   
   0-------3
   |       |
   |       |
   |       |
   1-------2 
*/
real MarchingCubes::surface(const Sign sig, const Comp & mcomp,
                            const int i, const int j, const int k) {
  return (cellface_covered(sig,mcomp,i,j,k)).value;
}

Heaviside::VAL2D MarchingCubes::cellface_covered(const Sign sig, const Comp & mcomp,
                                                 const int i, const int j, const int k) {
  CELL2D face;
  int isum(0);
  real surf;
  XY fcenter;

  if(mcomp == Comp::i()) {
    surf = (*clr).dSx(sig,i,j,k);
    fcenter.x = clr->yc(j);
    fcenter.y = clr->zc(k);
#if 0
    int ii = i;
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

      /* x: y->x, z->y */
      face.p[m].x = (*clr).yn(jj+1);
      face.p[m].y = (*clr).zn(kk+1);
#else
    int ii(i+1);
    if(sig == Sign::neg())
      ii = i;
    for(int m = 0; m != 4; m++) {
      int jj, kk;
      switch(m) {
        case(0) : jj = j  ; kk = k+1; break;
        case(1) : jj = j  ; kk = k  ; break;
        case(2) : jj = j+1; kk = k  ; break;
        case(3) : jj = j+1; kk = k+1; break;
      }

      face.val[m] = nodalvals[ii][jj][kk];

      /* x: y->x, z->y */
      face.p[m].x = (*clr).yn(jj);
      face.p[m].y = (*clr).zn(kk);
#endif

      if(fabs(face.val[m]-clrsurf) < boil::nano) {
        face.val[m]=clrsurf+boil::nano;
      }
     
      if(face.val[m] > clrsurf)
        isum++;
    }

    /* face-centered value */
    face.refval = 0.5*((*clr)[i][j][k]+(*clr)[ii-1][j][k]);
  }

  else if(mcomp == Comp::j()) {
    surf = (*clr).dSy(sig,i,j,k);
    fcenter.x = clr->xc(i);
    fcenter.y = clr->zc(k);
#if 0
    int jj = j;
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

      /* y: x->x, z->y */
      face.p[m].x = (*clr).xn(ii+1);
      face.p[m].y = (*clr).zn(kk+1);
#else
    int jj(j+1);
    if(sig == Sign::neg())
      jj = j;
    for(int m = 0; m != 4; m++) {
      int ii, kk;
      switch(m) {
        case(0) : ii = i  ; kk = k+1; break;
        case(1) : ii = i  ; kk = k  ; break;
        case(2) : ii = i+1; kk = k  ; break;
        case(3) : ii = i+1; kk = k+1; break;
      }

      face.val[m] = nodalvals[ii][jj][kk];

      /* y: x->x, z->y */
      face.p[m].x = (*clr).xn(ii);
      face.p[m].y = (*clr).zn(kk);
#endif

      if(fabs(face.val[m]-clrsurf) < boil::nano) {
        face.val[m]=clrsurf+boil::nano; 
      }

      if(face.val[m] > clrsurf)
        isum++;
    }

    /* face-centered value */
    face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i][jj-1][k]);
  }

  else if(mcomp == Comp::k()) {
    surf = (*clr).dSz(sig,i,j,k);
    fcenter.x = clr->xc(i);
    fcenter.y = clr->yc(j);
#if 0
    int kk = k;
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

      /* z: x->x, y->y */
      face.p[m].x = (*clr).xn(ii+1);
      face.p[m].y = (*clr).yn(jj+1);
#else
    int kk(k+1);
    if(sig == Sign::neg()) { 
      kk = k;
    }
    for(int m = 0; m != 4; m++) {
      int ii, jj;
      switch(m) {
        case(0) : ii = i  ; jj = j+1; break;
        case(1) : ii = i  ; jj = j  ; break;
        case(2) : ii = i+1; jj = j  ; break;
        case(3) : ii = i+1; jj = j+1; break;
      }

      face.val[m] = nodalvals[ii][jj][kk];

      /* z: x->x, y->y */
      face.p[m].x = (*clr).xn(ii);
      face.p[m].y = (*clr).yn(jj);
#endif

      if(fabs(face.val[m]-clrsurf) < boil::nano) {
        face.val[m]=clrsurf+boil::nano; 
      }

      if(face.val[m] > clrsurf)
        isum++;
    }

    /* face-centered value */
    face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i][j][kk-1]);
  }
 
  real af(0.0);
  if       (isum == 4) {
    af = 1.0;
  } else if(isum > 0) {
    std::vector<LINE> lines; /* dummy */
    XY centroid; /* dummy */
    af = 1.0-standing_square(face,clrsurf,surf,fcenter,lines,centroid)/surf;
  }

  VAL2D faceval;
  faceval.cell = face;
  faceval.value = af;
  return faceval;
}
