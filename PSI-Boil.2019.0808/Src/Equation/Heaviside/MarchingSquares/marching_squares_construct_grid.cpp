#include "marching_squares.h"

/******************************************************************************/
int MarchingSquares::construct_grid(const int i, const int j, const int k,
                                    CELL2D & grid) {
/***************************************************************************//**
*  \brief Fill the 2D grid for a given cell
*******************************************************************************/

/* 
   0-------3
   |       |
   |       |
   |       |
   1-------2
*/

  grid.refval = (*clr)[i][j][k];

  if(perpendicular == Comp::i()) {
    for(int m=0; m<=3; m++) {
      int ii,jj,kk;
      ii = i;
      switch(m) {
        case(0) : jj = j  ; kk = k+1; break;
        case(1) : jj = j  ; kk = k  ; break;
        case(2) : jj = j+1; kk = k  ; break;
        case(3) : jj = j+1; kk = k+1; break;
      }

      grid.val[m] = nodalvals[ii][jj][kk]; 
      grid.p[m].x = (*clr).yn(jj);
      grid.p[m].y = (*clr).zn(kk);
    }
  }

  else if(perpendicular == Comp::j()) {
    for(int m=0; m<=3; m++) {
      int ii,jj,kk;
      jj = j;
      switch(m) { 
        case(0) : ii = i  ; kk = k+1; break;
        case(1) : ii = i  ; kk = k  ; break;
        case(2) : ii = i+1; kk = k  ; break;
        case(3) : ii = i+1; kk = k+1; break;
      }

      grid.val[m] = nodalvals[ii][jj][kk];
      grid.p[m].x = (*clr).xn(ii);
      grid.p[m].y = (*clr).zn(kk);
    }
  }

  else if(perpendicular == Comp::k()) {
    for(int m=0; m<=3; m++) {
      int ii,jj,kk;
      kk = k;
      switch(m) {
        case(0) : ii = i  ; jj = j+1; break;
        case(1) : ii = i  ; jj = j  ; break;
        case(2) : ii = i+1; jj = j  ; break;
        case(3) : ii = i+1; jj = j+1; break;
      }

      grid.val[m] = nodalvals[ii][jj][kk];
      grid.p[m].x = (*clr).xn(ii);
      grid.p[m].y = (*clr).yn(jj);
    }
  }
 
  else {
    boil::aout<<"Marching Squares direction not properly set! Exiting."
              <<boil::endl;
    exit(0);
  }

  int isum(0);
  for(int m=0; m<=3; m++) {
    if(fabs(grid.val[m]-clrsurf) < boil::nano) {
#if 1
      if       ((*clr)[i][j][k]>(1.0-boil::nano)) {
        grid.val[m]=clrsurf-boil::nano;
      } else if((*clr)[i][j][k]<boil::nano) {
        grid.val[m]=clrsurf+boil::nano;
      } else {
        grid.val[m]=clrsurf+boil::nano;
      }
#else
      grid.val[m]=clrsurf-boil::nano;
#endif
    }

    if(grid.val[m]>clrsurf) isum++;
  }

#if 0
  boil::oout<<i<<" "<<j<<" "<<k<<" | ";
  for(int m(0); m<4;++m)
    boil::oout<<grid.val[m]<<" ";
  boil::oout<<boil::endl;
#endif

  return isum;
}
