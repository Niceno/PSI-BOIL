#include "marching_cubes.h"
#include <iostream>

/* calculate the surface area below a given isosurface */
real MarchingCubes::surfval(CELL2D grid, real isolevel) {
  int squareindex(0);
  real area(0.0);
  
  if (grid.val[0] < isolevel) squareindex |= 1; /* 0001 */
  if (grid.val[1] < isolevel) squareindex |= 2; /* 0010 */
  if (grid.val[2] < isolevel) squareindex |= 4; /* 0100 */
  if (grid.val[3] < isolevel) squareindex |= 8; /* 1000 */

  /* 0000: Face is entirely above the surface */
  
  /* 1111: Face is entirely below the surface */
  if (squareindex == 15)
    area += SurfaceArea4(grid.p[0],grid.p[1],grid.p[2],grid.p[3]);

  /* only one vertex is below the surface */

  /* 0001 */
  else if (squareindex == 1) {
    XY vertlist[3];
    
    vertlist[0] = grid.p[0];
    vertlist[1] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    vertlist[2] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

    area += SurfaceArea3(vertlist[0],vertlist[1],vertlist[2]);
  }
  /* 0010 */
  else if (squareindex == 2) {
    XY vertlist[3];

    vertlist[0] = grid.p[1];
    vertlist[1] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);
    vertlist[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);

    area += SurfaceArea3(vertlist[0],vertlist[1],vertlist[2]);
  }
  /* 0100 */
  else if (squareindex == 4) {
    XY vertlist[3];

    vertlist[0] = grid.p[2];
    vertlist[1] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);
    vertlist[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);

    area += SurfaceArea3(vertlist[0],vertlist[1],vertlist[2]);
  }
  /* 1000 */
  else if (squareindex == 8) {
    XY vertlist[3];

    vertlist[0] = grid.p[3];
    vertlist[1] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

    area += SurfaceArea3(vertlist[0],vertlist[1],vertlist[2]);
  }

  /* two adjacent vertices are below the surface */

  /* 0011 */
  else if (squareindex == 3) {
    XY vertlist[4];

    vertlist[0] = grid.p[0];
    vertlist[1] = grid.p[1];
    vertlist[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
    vertlist[3] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

    area += SurfaceArea4(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
  }
  /* 0110 */
  else if (squareindex == 6) {
    XY vertlist[4];

    vertlist[0] = grid.p[1];
    vertlist[1] = grid.p[2];
    vertlist[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
    vertlist[3] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);

    area += SurfaceArea4(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
  }
  /* 1001 */
  else if (squareindex == 9) {
    XY vertlist[4];

    vertlist[0] = grid.p[0];
    vertlist[1] = grid.p[3];
    vertlist[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);
    vertlist[3] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);

    area += SurfaceArea4(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
  }
  /* 1100 */
  else if (squareindex == 12) {
    XY vertlist[4];

    vertlist[0] = grid.p[2];
    vertlist[1] = grid.p[3];
    vertlist[2] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist[3] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);

    area += SurfaceArea4(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
  }

  /* two non-adjacent vertices are below the surface */

  /* 0101 */
  else if (squareindex == 5) {
    XY vertlist1[3];
    XY vertlist2[3];

    vertlist1[0] = grid.p[0];
    vertlist1[1] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    vertlist1[2] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

    area += SurfaceArea3(vertlist1[0],vertlist1[1],vertlist1[2]);

    vertlist2[0] = grid.p[2];
    vertlist2[1] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);
    vertlist2[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);

    area += SurfaceArea3(vertlist2[0],vertlist2[1],vertlist2[2]);
  }
  /* 1010 */
  else if (squareindex == 10) {
    XY vertlist1[3];
    XY vertlist2[3];

    vertlist1[0] = grid.p[1];
    vertlist1[1] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);
    vertlist1[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);

    area += SurfaceArea3(vertlist1[0],vertlist1[1],vertlist1[2]);

    vertlist2[0] = grid.p[3];
    vertlist2[1] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist2[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

    area += SurfaceArea3(vertlist2[0],vertlist2[1],vertlist2[2]);
  }

  /* three vertices are below the surface */

  /* 0111 */
  else if (squareindex == 7) {
    XY vertlist[5];

    vertlist[0] = grid.p[0];
    vertlist[1] = grid.p[1];
    vertlist[2] = grid.p[2];
    vertlist[3] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
    vertlist[4] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

    area += SurfaceArea5(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
  }
  /* 1011 */
  else if (squareindex == 11) {
    XY vertlist[5];

    vertlist[0] = grid.p[3];
    vertlist[1] = grid.p[0];
    vertlist[2] = grid.p[1];
    vertlist[3] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
    vertlist[4] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

    area += SurfaceArea5(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
  }
  /* 1101 */
  else if (squareindex == 13) {
    XY vertlist[5];

    vertlist[0] = grid.p[2];
    vertlist[1] = grid.p[3];
    vertlist[2] = grid.p[0];
    vertlist[3] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    vertlist[4] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);

    area += SurfaceArea5(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
  }
  /* 1110 */
  else if (squareindex == 14) {
    XY vertlist[5];

    vertlist[0] = grid.p[1];
    vertlist[1] = grid.p[2];
    vertlist[2] = grid.p[3];
    vertlist[3] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist[4] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);

    area += SurfaceArea5(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
  }

  return(area);
}

