#include "marching_squares.h"

/******************************************************************************/
real MarchingSquares::vf(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief calculate fraction of area of cell (i,j,k) above isoline
*******************************************************************************/

  if(dom->ibody().off(i,j,k))
    return 0.0;

  CELL2D grid;
  int isum = construct_grid(i,j,k,grid);

  real af(0.0);
  if       (isum == 4) {
    af = 1.0;
  } else if(isum > 0) {
    /* since we are working under the assumption of 2D cartesian, the
       following implementation is the only one safe */
    real surf;
    if       (perpendicular==Comp::i()) {
      surf = clr->dyc(j)*clr->dzc(k);
    } else if(perpendicular==Comp::j()) {
      surf = clr->dxc(i)*clr->dzc(k);
    } else if(perpendicular==Comp::k()) {
      surf = clr->dxc(i)*clr->dyc(j);
    } else {
      boil::aout<<"Marching Squares direction not properly set! Exiting."
                <<boil::endl;
      exit(0);
    }

    std::vector<LINE> lines; /* dummy */
    XY centroid;
    real are = standing_square(grid,clrsurf,surf,
                               {clr->xc(i),clr->zc(k)},
                               lines,centroid);
    /* whole function above evaluates area *below* IS! */
    af = 1.0-ratio(are,surf,centroid.x,clr->xc(i));
  }

  return af;
}
