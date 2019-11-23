#include "marching_cubes.h"

/******************************************************************************/
real MarchingCubes::vf(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief calculate fraction of volume of cell (i,j,k) above isosurface
*******************************************************************************/

  if(dom->ibody().off(i,j,k))
    return 0.0;

  CELL3D grid;
  int isum = construct_grid(i,j,k,grid);

  if(isum==8)
    return 1.0;
  if(isum==0)
    return 0.0;

  /* to achieve symmetry, cases isum < 4 are solved using an inverse problem */
  bool swtch(false);
  if(isum < 4) {
    swtch = true;
    for(int m=0; m<=7; m++) {
      grid.val[m] = 1.0 - grid.val[m];	
    }
  }

  real vol(0.0);
  vol += polygonise_volume(grid, clrsurf);

  std::vector<LINE> lines; /* dummy */
  for(int m = 0; m != 6; m++) {
    /* x: y->x, z->y */
    if(m==0) { /* west */
      CELL2D face;
      face.val[0] = grid.val[4]; face.p[0].x = grid.p[4].y; face.p[0].y = grid.p[4].z;
      face.val[1] = grid.val[0]; face.p[1].x = grid.p[0].y; face.p[1].y = grid.p[0].z;
      face.val[2] = grid.val[3]; face.p[2].x = grid.p[3].y; face.p[2].y = grid.p[3].z;
      face.val[3] = grid.val[7]; face.p[3].x = grid.p[7].y; face.p[3].y = grid.p[7].z;

      /* face-centered value */
      if(swtch) {
        face.refval = 1.0-0.5*((*clr)[i][j][k]+(*clr)[i-1][j][k]);
      } else {
        face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i-1][j][k]);
      }
     
      vol += (-1.0) * grid.p[4].x 
             * standing_square(face, clrsurf,clr->dSx(Sign::neg(),i,j,k),lines);
    }
    else if(m==1) { /* east */
      CELL2D face;
      face.val[0] = grid.val[5]; face.p[0].x = grid.p[5].y; face.p[0].y = grid.p[5].z;
      face.val[1] = grid.val[1]; face.p[1].x = grid.p[1].y; face.p[1].y = grid.p[1].z;
      face.val[2] = grid.val[2]; face.p[2].x = grid.p[2].y; face.p[2].y = grid.p[2].z;
      face.val[3] = grid.val[6]; face.p[3].x = grid.p[6].y; face.p[3].y = grid.p[6].z;

      /* face-centered value */
      if(swtch) {
        face.refval = 1.0-0.5*((*clr)[i][j][k]+(*clr)[i+1][j][k]);
      } else {
        face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i+1][j][k]);
      }
     
      vol += (+1.0) * grid.p[5].x
             * standing_square(face, clrsurf,clr->dSx(Sign::pos(),i,j,k),lines);
    }
    /* y: x->x, z->y */
    else if(m==2) { /* south */
      CELL2D face;
      face.val[0] = grid.val[4]; face.p[0].x = grid.p[4].x; face.p[0].y = grid.p[4].z;
      face.val[1] = grid.val[0]; face.p[1].x = grid.p[0].x; face.p[1].y = grid.p[0].z;
      face.val[2] = grid.val[1]; face.p[2].x = grid.p[1].x; face.p[2].y = grid.p[1].z;
      face.val[3] = grid.val[5]; face.p[3].x = grid.p[5].x; face.p[3].y = grid.p[5].z;

      /* face-centered value */
      if(swtch) {
        face.refval = 1.0-0.5*((*clr)[i][j][k]+(*clr)[i][j-1][k]);
      } else {
        face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i][j-1][k]);
      }
    
      vol += (-1.0) * grid.p[4].y
             * standing_square(face, clrsurf,clr->dSy(Sign::neg(),i,j,k),lines);
    }
    else if(m==3) { /* north */
      CELL2D face;
      face.val[0] = grid.val[7]; face.p[0].x = grid.p[7].x; face.p[0].y = grid.p[7].z;
      face.val[1] = grid.val[3]; face.p[1].x = grid.p[3].x; face.p[1].y = grid.p[3].z;
      face.val[2] = grid.val[2]; face.p[2].x = grid.p[2].x; face.p[2].y = grid.p[2].z;
      face.val[3] = grid.val[6]; face.p[3].x = grid.p[6].x; face.p[3].y = grid.p[6].z;

      /* face-centered value */
      if(swtch) {
        face.refval = 1.0-0.5*((*clr)[i][j][k]+(*clr)[i][j+1][k]);
      } else {
        face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i][j+1][k]);
      }
     
      vol += (+1.0) * grid.p[7].y
             * standing_square(face, clrsurf,clr->dSy(Sign::pos(),i,j,k),lines);
    }
    /* z: x->x, y->y */
    else if(m==4) { /* bottom */
      CELL2D face;
      face.val[0] = grid.val[3]; face.p[0].x = grid.p[3].x; face.p[0].y = grid.p[3].y;
      face.val[1] = grid.val[0]; face.p[1].x = grid.p[0].x; face.p[1].y = grid.p[0].y;
      face.val[2] = grid.val[1]; face.p[2].x = grid.p[1].x; face.p[2].y = grid.p[1].y;
      face.val[3] = grid.val[2]; face.p[3].x = grid.p[2].x; face.p[3].y = grid.p[2].y;

      /* face-centered value */
      if(swtch) {
        face.refval = 1.0-0.5*((*clr)[i][j][k]+(*clr)[i][j][k-1]);
      } else {
        face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i][j][k-1]);
      }
     
      vol += (-1.0) * grid.p[3].z
             * standing_square(face, clrsurf,clr->dSz(Sign::neg(),i,j,k),lines);
    }
    else if(m==5) { /* top */
      CELL2D face;
      face.val[0] = grid.val[7]; face.p[0].x = grid.p[7].x; face.p[0].y = grid.p[7].y;
      face.val[1] = grid.val[4]; face.p[1].x = grid.p[4].x; face.p[1].y = grid.p[4].y;
      face.val[2] = grid.val[5]; face.p[2].x = grid.p[5].x; face.p[2].y = grid.p[5].y;
      face.val[3] = grid.val[6]; face.p[3].x = grid.p[6].x; face.p[3].y = grid.p[6].y;

      /* face-centered value */
      if(swtch) {
        face.refval = 1.0-0.5*((*clr)[i][j][k]+(*clr)[i][j][k+1]);
      } else {
        face.refval = 0.5*((*clr)[i][j][k]+(*clr)[i][j][k+1]);
      }
     
      vol += (+1.0) * grid.p[7].z
             * standing_square(face, clrsurf,clr->dSy(Sign::pos(),i,j,k),lines);
    }
  }
  vol = vol/3. /clr->dV(i,j,k);

  if(vol > 1.0) {
    boil::aout << "marching_cubes::volume warning! color out of bounds at: "
               << i << " " << j << " " << k << " " << vol <<"\n";
    vol = 1.0;
  }
  if(vol < 0.0) {
    boil::aout << "marching_cubes::volume warning! color out of bounds at: "
               << i << " " << j << " " << k << " " << vol <<"\n";
    vol = 0.0;
  }

  /* whole function above evaluates volume *below* IS; flagging is inverted! */
  if(swtch)
    return vol;  
  else
    return (1.0-vol);
}
