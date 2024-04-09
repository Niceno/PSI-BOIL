#include "marching_cube.h"

/******************************************************************************/
real MarchingCube::volume(const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate fraction of volume of cell (i,j,k) below isosurface
*******************************************************************************/

  if(dom->ibody().off(i,j,k)) return 0.0;

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
    grid.val[m]=((*clr)[ii][jj  ][kk  ]+(*clr)[ii+1][jj  ][kk  ]
                +(*clr)[ii][jj+1][kk  ]+(*clr)[ii+1][jj+1][kk  ]
                +(*clr)[ii][jj  ][kk+1]+(*clr)[ii+1][jj  ][kk+1]
                +(*clr)[ii][jj+1][kk+1]+(*clr)[ii+1][jj+1][kk+1])/8.0;
#if 1
    if(fabs(grid.val[m] - clrsurf) < boil::pico){
        grid.val[m]=clrsurf+boil::pico;
    }
#else
    /* for the purpose of indicator reconstruction, the sign of phase change
     * should not matter? Cf. phasechange_marching_cube.cpp */
#endif
    if(grid.val[m]>clrsurf)isum++;
  }
  if(isum==8)return(1.0);
  if(isum==0)return(0.0);
  
  bool swtch(false);
  
  /* to achieve symmetry, cases isum < 4 are solved using an inverse problem */
  if (isum < 4) {
    swtch = true;
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

  real volume(0.0);
  volume += polygonise_volume(grid, clrsurf);
  for (int m = 0; m != 6; m++) {
    /* x: y->x, z->y */
    if(m==0) { /* west */
      CELLFACE face;
      face.val[0] = grid.val[4]; face.p[0].x = grid.p[4].y; face.p[0].y = grid.p[4].z;
      face.val[1] = grid.val[0]; face.p[1].x = grid.p[0].y; face.p[1].y = grid.p[0].z;
      face.val[2] = grid.val[3]; face.p[2].x = grid.p[3].y; face.p[2].y = grid.p[3].z;
      face.val[3] = grid.val[7]; face.p[3].x = grid.p[7].y; face.p[3].y = grid.p[7].z;
     
      volume += (-1.0) * grid.p[4].x * surfval(face, clrsurf);
    }
    else if(m==1) { /* east */
      CELLFACE face;
      face.val[0] = grid.val[5]; face.p[0].x = grid.p[5].y; face.p[0].y = grid.p[5].z;
      face.val[1] = grid.val[1]; face.p[1].x = grid.p[1].y; face.p[1].y = grid.p[1].z;
      face.val[2] = grid.val[2]; face.p[2].x = grid.p[2].y; face.p[2].y = grid.p[2].z;
      face.val[3] = grid.val[6]; face.p[3].x = grid.p[6].y; face.p[3].y = grid.p[6].z;
     
      volume += (+1.0) * grid.p[5].x * surfval(face, clrsurf);
    }
    /* y: x->x, z->y */
    else if(m==2) { /* south */
      CELLFACE face;
      face.val[0] = grid.val[4]; face.p[0].x = grid.p[4].x; face.p[0].y = grid.p[4].z;
      face.val[1] = grid.val[0]; face.p[1].x = grid.p[0].x; face.p[1].y = grid.p[0].z;
      face.val[2] = grid.val[1]; face.p[2].x = grid.p[1].x; face.p[2].y = grid.p[1].z;
      face.val[3] = grid.val[5]; face.p[3].x = grid.p[5].x; face.p[3].y = grid.p[5].z;
    
      volume += (-1.0) * grid.p[4].y * surfval(face, clrsurf);
    }
    else if(m==3) { /* north */
      CELLFACE face;
      face.val[0] = grid.val[7]; face.p[0].x = grid.p[7].x; face.p[0].y = grid.p[7].z;
      face.val[1] = grid.val[3]; face.p[1].x = grid.p[3].x; face.p[1].y = grid.p[3].z;
      face.val[2] = grid.val[2]; face.p[2].x = grid.p[2].x; face.p[2].y = grid.p[2].z;
      face.val[3] = grid.val[6]; face.p[3].x = grid.p[6].x; face.p[3].y = grid.p[6].z;
     
      volume += (+1.0) * grid.p[7].y * surfval(face, clrsurf);
    }
    /* z: x->x, y->y */
    else if(m==4) { /* bottom */
      CELLFACE face;
      face.val[0] = grid.val[3]; face.p[0].x = grid.p[3].x; face.p[0].y = grid.p[3].y;
      face.val[1] = grid.val[0]; face.p[1].x = grid.p[0].x; face.p[1].y = grid.p[0].y;
      face.val[2] = grid.val[1]; face.p[2].x = grid.p[1].x; face.p[2].y = grid.p[1].y;
      face.val[3] = grid.val[2]; face.p[3].x = grid.p[2].x; face.p[3].y = grid.p[2].y;
     
      volume += (-1.0) * grid.p[3].z * surfval(face, clrsurf);
    }
    else if(m==5) { /* top */
      CELLFACE face;
      face.val[0] = grid.val[7]; face.p[0].x = grid.p[7].x; face.p[0].y = grid.p[7].y;
      face.val[1] = grid.val[4]; face.p[1].x = grid.p[4].x; face.p[1].y = grid.p[4].y;
      face.val[2] = grid.val[5]; face.p[2].x = grid.p[5].x; face.p[2].y = grid.p[5].y;
      face.val[3] = grid.val[6]; face.p[3].x = grid.p[6].x; face.p[3].y = grid.p[6].y;
     
      volume += (+1.0) * grid.p[7].z * surfval(face, clrsurf);
    }
  }
  volume = volume/3 /(*clr).dV(i,j,k);

  if (volume > 1.0) {
    boil::oout << "marching_cube::volume warning! color out of bounds at: "
               << i << " " << j << " " << k << " " << volume <<"\n";
    volume = 1.0;
  }
  if (volume < 0.0) {
    boil::oout << "marching_cube::volume warning! color out of bounds at: "
               << i << " " << j << " " << k << " " << volume <<"\n";
    volume = 0.0;
  }

  if (swtch) return (volume);  
  else return (1.0-volume);
}
