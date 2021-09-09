#include "marching_cubes.h"

/******************************************************************************/
real MarchingCubes::ad(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief calculate iso-surface area density of cell (i,j,k)
*******************************************************************************/

  if(dom->ibody().off(i,j,k))
    return 0.0;

  CELL3D grid;
  int isum = construct_grid(i,j,k,grid);

  if(isum==0||isum==8)
    return 0.0;

  /* to achieve symmetry, cases isum < 4 are solved using an inverse problem */
  if(isum < 4) {
    for(int m=0; m<=7; m++) {
      grid.val[m] = 1.0 - grid.val[m];
    }
  }

  real areaval = polygonise_area(grid, clrsurf);

  return areaval/clr->dV(i,j,k);
}
