#include "marching_squares.h"
  
/******************************************************************************/
int MarchingSquares::status(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief check status of cell (i,j,k)
          -1: fully below the isosurface
           0: at the interface
          +1: fully above the isosurface
       -1000: immersed body
*******************************************************************************/

  if(dom->ibody().off(i,j,k))
    return -1000;

  CELL2D grid;
  int isum = construct_grid(i,j,k,grid);

  if      (isum==0) {
    return -1;
  } else if(isum==4) {
    return +1;
  } else {
    return 0;
  }
}
