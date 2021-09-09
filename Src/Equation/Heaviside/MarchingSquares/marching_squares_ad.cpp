#include "marching_squares.h"

/******************************************************************************/
real MarchingSquares::ad(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief calculate iso-line length density of cell (i,j,k)
*******************************************************************************/
  std::vector<LINE> lines; /* dummy */
  return ad(i,j,k,lines);
}

/******************************************************************************/
real MarchingSquares::ad(const int i, const int j, const int k,
                         std::vector<LINE> & lines) {
/******************************************************************************/
  lines.clear();

  if(dom->ibody().off(i,j,k))
    return 0.0;

  CELL2D grid;
  int isum = construct_grid(i,j,k,grid);

  if(isum>0&&isum<4) {
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

    standing_square(grid,clrsurf,surf,lines);

    /* allow for virtual implementation */
    return line_density(lines,surf,clr->xc(i));
  }

  return 0.0;
}
