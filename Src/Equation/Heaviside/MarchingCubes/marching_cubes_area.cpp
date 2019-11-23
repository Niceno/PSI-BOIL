#include "marching_cubes.h"

/******************************************************************************/
real MarchingCubes::area(const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate iso-surface area of cell (i,j,k)
*******************************************************************************/

  if(dom->ibody().off(i,j,k))
    return 0.0;

  CELL3D grid;

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

    if(fabs(grid.val[m]-clrsurf) < boil::pico) {
#if 1
      if       ((*clr)[i][j][k]>(1.0-boil::pico)) {
        grid.val[m]=clrsurf-boil::pico;
      } else if((*clr)[i][j][k]<boil::pico) {
        grid.val[m]=clrsurf+boil::pico;
      } else {
        grid.val[m]=clrsurf+boil::pico;
      }
#else
      grid.val[m]=clrsurf-boil::pico;
#endif
    }

    if(grid.val[m]>clrsurf) isum++;
  }

  if(isum==0||isum==8)
    return(0.0);

  /* to achieve symmetry, cases isum < 4 are solved using an inverse problem */
  if (isum < 4) {
    for (int m=0; m<=7; m++) {
      grid.val[m] = 1.0 - grid.val[m];
    }
  }

  real areaval = polygonise_area(grid, clrsurf);

#if 0
  if(i==3&&j==3&&k==7) {
    for(int i(0); i < 8; ++i)
      boil::oout<<grid.val[i]<<boil::endl;
    boil::oout<<areaval<<boil::endl;
  }
#endif

  return (areaval);
}
