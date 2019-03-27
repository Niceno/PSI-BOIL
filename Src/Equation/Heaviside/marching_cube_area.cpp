#include "marching_cube.h"

/******************************************************************************/
real MarchingCube::area(const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate iso-surface area of cell (i,j,k)
*******************************************************************************/
  GRIDCELL grid;

  int isum=0;
  for (int m=0; m<=7; m++){
    int ii,jj,kk;
    if(m==0)     {ii=i-1; jj=j-1; kk=k-1;}
    else if(m==1){ii=i  ; jj=j-1; kk=k-1;}
    else if(m==2){ii=i  ; jj=j  ; kk=k-1;}
    else if(m==3){ii=i-1; jj=j  ; kk=k-1;}
    else if(m==4){ii=i-1; jj=j-1; kk=k  ;}
    else if(m==5){ii=i  ; jj=j-1; kk=k  ;}
    else if(m==6){ii=i  ; jj=j  ; kk=k  ;}
    else if(m==7){ii=i-1; jj=j  ; kk=k  ;}
#if 0
    grid.val[m]=((*clr)[ii][jj  ][kk  ]+(*clr)[ii+1][jj  ][kk  ]
                +(*clr)[ii][jj+1][kk  ]+(*clr)[ii+1][jj+1][kk  ]
                +(*clr)[ii][jj  ][kk+1]+(*clr)[ii+1][jj  ][kk+1]
                +(*clr)[ii][jj+1][kk+1]+(*clr)[ii+1][jj+1][kk+1])/8.0;
#else
    grid.val[m]=0.0;
    for(int idx=0;idx<2;idx++)
      for(int jdx=0;jdx<2;jdx++)
        for(int kdx=0;kdx<2;kdx++)
          grid.val[m]+=std::max(0.0,std::min(1.0,(*clr)[ii+idx][jj+jdx][kk+kdx]));
    grid.val[m] /= 8.0;
#endif 


    if(fabs(grid.val[m] - clrsurf) < boil::pico){
#if 1
      if( (*clr)[i][j][k]>(1.0-boil::pico)){
        grid.val[m]=clrsurf-boil::pico;
      } else if ( (*clr)[i][j][k]<(boil::pico)){
        grid.val[m]=clrsurf+boil::pico;
      } else {
        grid.val[m]=clrsurf+boil::pico;
      }
#else
      grid.val[m]=clrsurf-boil::pico;
#endif
    }

    if(grid.val[m]>clrsurf)isum++;
  }

  if(isum==0||isum==8)return(0.0);

  /* to achieve symmetry, cases isum < 4 are solved using an inverse problem */
  if (isum < 4) {
    for (int m=0; m<=7; m++) {
      grid.val[m] = 1.0 - grid.val[m];
    }
  }

  for (int m=0; m<=7; m++){
    int ii,jj,kk;
    if(m==0)     {ii=i  ; jj=j  ; kk=k  ;}
    else if(m==1){ii=i+1; jj=j  ; kk=k  ;}
    else if(m==2){ii=i+1; jj=j+1; kk=k  ;}
    else if(m==3){ii=i  ; jj=j+1; kk=k  ;}
    else if(m==4){ii=i  ; jj=j  ; kk=k+1;}
    else if(m==5){ii=i+1; jj=j  ; kk=k+1;}
    else if(m==6){ii=i+1; jj=j+1; kk=k+1;}
    else if(m==7){ii=i  ; jj=j+1; kk=k+1;}
    grid.p[m].x=(*clr).xn(ii);
    grid.p[m].y=(*clr).yn(jj);
    grid.p[m].z=(*clr).zn(kk);
  }

#if 0
  grid.p[0].x=0.0; grid.p[0].y=0.0; grid.p[0].z=0.0; grid.val[0]=0.5;
  grid.p[1].x=1.0; grid.p[1].y=0.0; grid.p[1].z=0.0; grid.val[1]=0.5;
  grid.p[2].x=1.0; grid.p[2].y=1.0; grid.p[2].z=0.0; grid.val[2]=0.5;
  grid.p[3].x=0.0; grid.p[3].y=1.0; grid.p[3].z=0.0; grid.val[3]=0.5;
  grid.p[4].x=0.0; grid.p[4].y=0.0; grid.p[4].z=1.0; grid.val[4]=0.0;
  grid.p[5].x=1.0; grid.p[5].y=0.0; grid.p[5].z=1.0; grid.val[5]=0.0;
  grid.p[6].x=1.0; grid.p[6].y=1.0; grid.p[6].z=1.0; grid.val[6]=0.0;
  grid.p[7].x=0.0; grid.p[7].y=1.0; grid.p[7].z=1.0; grid.val[7]=0.0;
#endif

  real area = polygonise_area(grid, clrsurf);
  //std::cout<<"marching_cube "<<area<<"\n";
  //exit(0);
 
  return (area);
}
