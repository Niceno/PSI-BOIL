#include "heaviside.h"

/******************************************************************************/
real Heaviside::standing_square(const CELL2D & grid, const real & isolevel,
                                const real & totarea,
                                std::vector<LINE> & lines) {
/***************************************************************************//**
*  \brief Core of the Marching Squares algorithm. Provides surface area under a
          given isosurface, together with line(s) defining the isosurface. This 
          function is used by both Marching Cubes and Marching Squares.
*******************************************************************************/

/* 
   0-------3
   |       |
   |       |
   |       |
   1-------2
*/

  lines.clear();

  int squareindex(0);
  real are(0.0);
  
  if (grid.val[0] < isolevel) squareindex |= 1; /* 0001 */
  if (grid.val[1] < isolevel) squareindex |= 2; /* 0010 */
  if (grid.val[2] < isolevel) squareindex |= 4; /* 0100 */
  if (grid.val[3] < isolevel) squareindex |= 8; /* 1000 */

  /* 0000: Face is entirely above the surface */
  
  /* 1111: Face is entirely below the surface */
  if (squareindex == 15)
    are = Shoelace(grid.p[0],grid.p[1],grid.p[2],grid.p[3]);

  /* only one vertex is below the surface */

  /* 0001 */
  else if (squareindex == 1) {
    XY vertlist[3];
    
    vertlist[0] = grid.p[0];
    vertlist[1] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    vertlist[2] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2]);
    lines.push_back({vertlist[1],vertlist[2]});
  }
  /* 0010 */
  else if (squareindex == 2) {
    XY vertlist[3];

    vertlist[0] = grid.p[1];
    vertlist[1] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);
    vertlist[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2]);
    lines.push_back({vertlist[1],vertlist[2]});
  }
  /* 0100 */
  else if (squareindex == 4) {
    XY vertlist[3];

    vertlist[0] = grid.p[2];
    vertlist[1] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);
    vertlist[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2]);
    lines.push_back({vertlist[1],vertlist[2]});
  }
  /* 1000 */
  else if (squareindex == 8) {
    XY vertlist[3];

    vertlist[0] = grid.p[3];
    vertlist[1] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2]);
    lines.push_back({vertlist[1],vertlist[2]});
  }

  /* two adjacent vertices are below the surface */

  /* 0011 */
  else if (squareindex == 3) {
    XY vertlist[4];

    vertlist[0] = grid.p[0];
    vertlist[1] = grid.p[1];
    vertlist[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
    vertlist[3] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
    lines.push_back({vertlist[2],vertlist[3]});
  }
  /* 0110 */
  else if (squareindex == 6) {
    XY vertlist[4];

    vertlist[0] = grid.p[1];
    vertlist[1] = grid.p[2];
    vertlist[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
    vertlist[3] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
    lines.push_back({vertlist[2],vertlist[3]});
  }
  /* 1001 */
  else if (squareindex == 9) {
    XY vertlist[4];

    vertlist[0] = grid.p[0];
    vertlist[1] = grid.p[3];
    vertlist[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);
    vertlist[3] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
    lines.push_back({vertlist[2],vertlist[3]});
  }
  /* 1100 */
  else if (squareindex == 12) {
    XY vertlist[4];

    vertlist[0] = grid.p[2];
    vertlist[1] = grid.p[3];
    vertlist[2] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist[3] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3]);
    lines.push_back({vertlist[2],vertlist[3]});
  }

  /* two non-adjacent vertices are below the surface */

  /* 0101 */
  else if (squareindex == 5) {
    XY vertlist1[3];
    XY vertlist2[3];

    /* we have to distinguish based on the cell centre */
    if(grid.refval>=isolevel) {
      /* northwest corner */
      vertlist1[0] = grid.p[0];
      /* between corner and its neighbor in 'south' direction */
      vertlist1[1] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
      /* between corner and its neighbor in 'east' direction */
      vertlist1[2] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

      are += Shoelace(vertlist1[0],vertlist1[1],vertlist1[2]);
      lines.push_back({vertlist1[1],vertlist1[2]});

      /* southeast corner */
      vertlist2[0] = grid.p[2];
      /* between corner and its neighbor in 'west' direction */
      vertlist2[1] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);
      /* between corner and its neighbor in 'north' direction */
      vertlist2[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);

      are += Shoelace(vertlist2[0],vertlist2[1],vertlist2[2]);
      lines.push_back({vertlist2[1],vertlist2[2]});
    } else {
      /* southwest corner */
      vertlist1[0] = grid.p[1];
      /* between corner and its neighbor in 'north' direction */
      vertlist1[1] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);
      /* between corner and its neighbor in 'east' direction */
      vertlist1[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);

      real complement = Shoelace(vertlist1[0],vertlist1[1],vertlist1[2]);
      lines.push_back({vertlist1[1],vertlist1[2]});
  
      /* northeast corner */
      vertlist2[0] = grid.p[3];
      /* between corner and its neighbor in 'west' direction */
      vertlist2[1] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
      /* between corner and its neighbor in 'south' direction */
      vertlist2[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

      complement += Shoelace(vertlist2[0],vertlist2[1],vertlist2[2]);
      lines.push_back({vertlist2[1],vertlist2[2]});
      
      /* we are solving the complementary problem (make a drawing...) */
      are = totarea - complement;
    }
  }
  /* 1010 */
  else if (squareindex == 10) {
    XY vertlist1[3];
    XY vertlist2[3];

    if(grid.refval>=isolevel) {
      vertlist1[0] = grid.p[1];
      vertlist1[1] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);
      vertlist1[2] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
  
      are += Shoelace(vertlist1[0],vertlist1[1],vertlist1[2]);
      lines.push_back({vertlist1[1],vertlist1[2]});
  
      vertlist2[0] = grid.p[3];
      vertlist2[1] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
      vertlist2[2] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

      are += Shoelace(vertlist2[0],vertlist2[1],vertlist2[2]);
      lines.push_back({vertlist2[1],vertlist2[2]});
    } else {
      vertlist1[0] = grid.p[0];
      vertlist1[1] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
      vertlist1[2] = VertexInterp(isolevel,grid.p[0],grid.p[3],grid.val[0],grid.val[3]);

      real complement = Shoelace(vertlist1[0],vertlist1[1],vertlist1[2]);
      lines.push_back({vertlist1[1],vertlist1[2]});

      vertlist2[0] = grid.p[2];
      vertlist2[1] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);
      vertlist2[2] = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);

      complement += Shoelace(vertlist2[0],vertlist2[1],vertlist2[2]);
      lines.push_back({vertlist2[1],vertlist2[2]});

      are = totarea - complement;
    }
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

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
    lines.push_back({vertlist[3],vertlist[4]});
  }
  /* 1011 */
  else if (squareindex == 11) {
    XY vertlist[5];

    vertlist[0] = grid.p[3];
    vertlist[1] = grid.p[0];
    vertlist[2] = grid.p[1];
    vertlist[3] = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
    vertlist[4] = VertexInterp(isolevel,grid.p[3],grid.p[2],grid.val[3],grid.val[2]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
    lines.push_back({vertlist[3],vertlist[4]});
  }
  /* 1101 */
  else if (squareindex == 13) {
    XY vertlist[5];

    vertlist[0] = grid.p[2];
    vertlist[1] = grid.p[3];
    vertlist[2] = grid.p[0];
    vertlist[3] = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    vertlist[4] = VertexInterp(isolevel,grid.p[2],grid.p[1],grid.val[2],grid.val[1]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
    lines.push_back({vertlist[3],vertlist[4]});
  }
  /* 1110 */
  else if (squareindex == 14) {
    XY vertlist[5];

    vertlist[0] = grid.p[1];
    vertlist[1] = grid.p[2];
    vertlist[2] = grid.p[3];
    vertlist[3] = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    vertlist[4] = VertexInterp(isolevel,grid.p[1],grid.p[0],grid.val[1],grid.val[0]);

    are = Shoelace(vertlist[0],vertlist[1],vertlist[2],vertlist[3],vertlist[4]);
    lines.push_back({vertlist[3],vertlist[4]});
  }

  return are;
}

