#include "marching_cubes.h"
  
/******************************************************************************/
int MarchingCubes::construct_grid(const int i, const int j, const int k,
                                  CELL3D & grid) {
/***************************************************************************//**
*  \brief Fill the 3D grid for a given cell
*******************************************************************************/
/*
   https://codegolf.stackexchange.com/questions/126644/draw-an-ascii-cuboid

      7--------6
     /|       /|
    / |      / |
   4--------5  |
   |  3-----|--2
   | /      | /
   |/       |/
   0--------1
*/

  int isum(0);
  for (int m=0; m<=7; m++) {
    int ii,jj,kk;
#if 0
    if(m==0)     {ii=i-1; jj=j-1; kk=k-1;}
    else if(m==1){ii=i  ; jj=j-1; kk=k-1;}
    else if(m==2){ii=i  ; jj=j  ; kk=k-1;}
    else if(m==3){ii=i-1; jj=j  ; kk=k-1;}
    else if(m==4){ii=i-1; jj=j-1; kk=k  ;}
    else if(m==5){ii=i  ; jj=j-1; kk=k  ;}
    else if(m==6){ii=i  ; jj=j  ; kk=k  ;}
    else if(m==7){ii=i-1; jj=j  ; kk=k  ;}

    grid.val[m]=0.0;
    for(int idx=0;idx<2;idx++)
      for(int jdx=0;jdx<2;jdx++)
        for(int kdx=0;kdx<2;kdx++)
          grid.val[m]+=std::max(0.0,std::min(1.0,(*clr)[ii+idx][jj+jdx][kk+kdx]));
    grid.val[m] /= 8.0;

    switch(m) {
      case(0) : ii = i  ; jj = j  ; kk = k  ; break;
      case(1) : ii = i+1; jj = j  ; kk = k  ; break;
      case(2) : ii = i+1; jj = j+1; kk = k  ; break;
      case(3) : ii = i  ; jj = j+1; kk = k  ; break;
      case(4) : ii = i  ; jj = j  ; kk = k+1; break;
      case(5) : ii = i+1; jj = j  ; kk = k+1; break;
      case(6) : ii = i+1; jj = j+1; kk = k+1; break;
      case(7) : ii = i  ; jj = j+1; kk = k+1; break;
    }
    grid.p[m].x = (*clr).xn(ii);
    grid.p[m].y = (*clr).yn(jj);
    grid.p[m].z = (*clr).zn(kk);
#else 
    switch(m) {
      case(0) : ii = i  ; jj = j  ; kk = k  ; break;
      case(1) : ii = i+1; jj = j  ; kk = k  ; break;
      case(2) : ii = i+1; jj = j+1; kk = k  ; break;
      case(3) : ii = i  ; jj = j+1; kk = k  ; break;
      case(4) : ii = i  ; jj = j  ; kk = k+1; break;
      case(5) : ii = i+1; jj = j  ; kk = k+1; break;
      case(6) : ii = i+1; jj = j+1; kk = k+1; break;
      case(7) : ii = i  ; jj = j+1; kk = k+1; break;
    }

    grid.val[m] = nodalvals[ii][jj][kk]; 
    grid.p[m].x = (*clr).xn(ii);
    grid.p[m].y = (*clr).yn(jj);
    grid.p[m].z = (*clr).zn(kk);
#endif 

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

  return isum;
}
