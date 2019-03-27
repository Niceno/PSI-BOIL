#include "marching_cube.h"
#include <iostream>
#include <iomanip>

/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
real MarchingCube::polygonise_volume(GRIDCELL grid, real isolevel)
{
   TRIANGLE triangles[5];
   int i,ntriang;
   int cubeindex;
   VERT vertlist[12];

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
   if (grid.val[0] < isolevel) cubeindex |= 1;
   if (grid.val[1] < isolevel) cubeindex |= 2;
   if (grid.val[2] < isolevel) cubeindex |= 4;
   if (grid.val[3] < isolevel) cubeindex |= 8;
   if (grid.val[4] < isolevel) cubeindex |= 16;
   if (grid.val[5] < isolevel) cubeindex |= 32;
   if (grid.val[6] < isolevel) cubeindex |= 64;
   if (grid.val[7] < isolevel) cubeindex |= 128;

   /* Cube is entirely above the surface */
   if (edgeTable_volume[cubeindex] == 0)
      return(1.0);

   /* Cube is entirely below the surface */
   if (edgeTable_volume[cubeindex] == 0x2000)
      return(0.0);

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable_volume[cubeindex] & 1) {
      vertlist[0].v =
         VertexInterp3D(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
      if (grid.val[0] < isolevel) vertlist[0].ref = grid.p[0];
      else vertlist[0].ref = grid.p[1];
   }
   if (edgeTable_volume[cubeindex] & 2) {
      vertlist[1].v =
         VertexInterp3D(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
      if (grid.val[1] < isolevel) vertlist[1].ref = grid.p[1];
      else vertlist[1].ref = grid.p[2];
   }
   if (edgeTable_volume[cubeindex] & 4) {
      vertlist[2].v =
         VertexInterp3D(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
      if (grid.val[2] < isolevel) vertlist[2].ref = grid.p[2];
      else vertlist[2].ref = grid.p[3];
   }
   if (edgeTable_volume[cubeindex] & 8) {
      vertlist[3].v =
         VertexInterp3D(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
      if (grid.val[3] < isolevel) vertlist[3].ref = grid.p[3];
      else vertlist[3].ref = grid.p[0];
   }
   if (edgeTable_volume[cubeindex] & 16) {
      vertlist[4].v =
         VertexInterp3D(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
      if (grid.val[4] < isolevel) vertlist[4].ref = grid.p[4];
      else vertlist[4].ref = grid.p[5];
   }
   if (edgeTable_volume[cubeindex] & 32) {
      vertlist[5].v =
         VertexInterp3D(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
      if (grid.val[5] < isolevel) vertlist[5].ref = grid.p[5];
      else vertlist[5].ref = grid.p[6];
   }
   if (edgeTable_volume[cubeindex] & 64) {
      vertlist[6].v =
         VertexInterp3D(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
      if (grid.val[6] < isolevel) vertlist[6].ref = grid.p[6];
      else vertlist[6].ref = grid.p[7];
   }
   if (edgeTable_volume[cubeindex] & 128) {
      vertlist[7].v =
         VertexInterp3D(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
      if (grid.val[7] < isolevel) vertlist[7].ref = grid.p[7];
      else vertlist[7].ref = grid.p[4];
   }
   if (edgeTable_volume[cubeindex] & 256) {
      vertlist[8].v =
         VertexInterp3D(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
      if (grid.val[0] < isolevel) vertlist[8].ref = grid.p[0];
      else vertlist[8].ref = grid.p[4];
   }
   if (edgeTable_volume[cubeindex] & 512) {
      vertlist[9].v =
         VertexInterp3D(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
      if (grid.val[1] < isolevel) vertlist[9].ref = grid.p[1];
      else vertlist[9].ref = grid.p[5];
   }
   if (edgeTable_volume[cubeindex] & 1024) {
      vertlist[10].v =
         VertexInterp3D(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
      if (grid.val[2] < isolevel) vertlist[10].ref = grid.p[2];
      else vertlist[10].ref = grid.p[6];
   }
   if (edgeTable_volume[cubeindex] & 2048) {
      vertlist[11].v =
         VertexInterp3D(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);
      if (grid.val[3] < isolevel) vertlist[11].ref = grid.p[3];
      else vertlist[11].ref = grid.p[7];
   }

   /* Create the triangle */
   ntriang = 0;
   for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
      triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]].v;
      triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]].v;
      triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]].v;

      triangles[ntriang].v[0] = vertlist[triTable[cubeindex][i  ]].ref;
      triangles[ntriang].v[1] = vertlist[triTable[cubeindex][i+1]].ref;
      triangles[ntriang].v[2] = vertlist[triTable[cubeindex][i+2]].ref;
      ntriang++;
   }

   real vol_area = 0.0;
   for (int itri=0; itri<ntriang; itri++){
      vol_area += triangle_vol_area(triangles[itri]);
   }

   return(vol_area);
}

/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
real MarchingCube::polygonise_area(GRIDCELL grid, real isolevel)
{
   TRIANGLE triangles[5];
   int i,ntriang;
   int cubeindex;
   XYZ vertlist[12];

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
   if (grid.val[0] < isolevel) cubeindex |= 1;
   if (grid.val[1] < isolevel) cubeindex |= 2;
   if (grid.val[2] < isolevel) cubeindex |= 4;
   if (grid.val[3] < isolevel) cubeindex |= 8;
   if (grid.val[4] < isolevel) cubeindex |= 16;
   if (grid.val[5] < isolevel) cubeindex |= 32;
   if (grid.val[6] < isolevel) cubeindex |= 64;
   if (grid.val[7] < isolevel) cubeindex |= 128;

   /* Cube is entirely in/out of the surface */
   if (edgeTable_area[cubeindex] == 0)
      return(0.0);

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable_area[cubeindex] & 1)
      vertlist[0] =
         VertexInterp3D(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
   if (edgeTable_area[cubeindex] & 2)
      vertlist[1] =
         VertexInterp3D(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
   if (edgeTable_area[cubeindex] & 4)
      vertlist[2] =
         VertexInterp3D(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
   if (edgeTable_area[cubeindex] & 8)
      vertlist[3] =
         VertexInterp3D(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
   if (edgeTable_area[cubeindex] & 16)
      vertlist[4] =
         VertexInterp3D(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
   if (edgeTable_area[cubeindex] & 32)
      vertlist[5] =
         VertexInterp3D(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
   if (edgeTable_area[cubeindex] & 64)
      vertlist[6] =
         VertexInterp3D(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
   if (edgeTable_area[cubeindex] & 128)
      vertlist[7] =
         VertexInterp3D(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
   if (edgeTable_area[cubeindex] & 256)
      vertlist[8] =
         VertexInterp3D(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
   if (edgeTable_area[cubeindex] & 512)
      vertlist[9] =
         VertexInterp3D(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
   if (edgeTable_area[cubeindex] & 1024)
      vertlist[10] =
         VertexInterp3D(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
   if (edgeTable_area[cubeindex] & 2048)
      vertlist[11] =
         VertexInterp3D(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

   /* Create the triangle */
   ntriang = 0;
   for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
      triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
      triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
      triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];
      ntriang++;
   }

   real area = 0.0;
   for (int itri=0; itri<ntriang; itri++){
      area += triangle_area(triangles[itri]);
   }

   return(area);
}
