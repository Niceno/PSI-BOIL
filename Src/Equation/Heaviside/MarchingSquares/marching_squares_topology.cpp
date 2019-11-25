#include "marching_squares.h"

/******************************************************************************/
void MarchingSquares::topology(Topology & topo) {
/***************************************************************************//**
*  \brief build the interface topology (except for iflag)
*******************************************************************************/

  /* step zero: update node values */
  evaluate_nodes();

  /* step one: adens, normal vector + line constant (in stmp) */
  std::vector<LINE> lines;

  /* since this is 2D, one normal component will always be zero. */
  Scalar * nnx;
  Scalar * nny;
  Scalar * nzero;

  if       (perpendicular==Comp::i()) {
    nnx = topo.ny;
    nny = topo.nz;
    nzero = topo.nx;
  } else if(perpendicular==Comp::j()) {
    nnx = topo.nx;
    nny = topo.nz;
    nzero = topo.ny;
  } else if(perpendicular==Comp::k()) {
    nnx = topo.nx;
    nny = topo.ny;
    nzero = topo.nz;
  } else {
    boil::aout<<"Marching Squares direction not properly set! Exiting."
              <<boil::endl;
    exit(0);
  }

  for_vijk(stmp,i,j,k) {
    /* lines get cleared with every call of this function */
    (*topo.adens) = ad(i,j,k,lines);

    extract_line_parameters(lines,(*nnx)[i][j][k],(*nny)[i][j][k],stmp[i][j][k]);
    (*nzero)[i][j][k] = 0.0;
  }

  nnx->exchange_all();
  nny->exchange_all();
  nzero->exchange_all();
  stmp.exchange_all();

  /* step two: calculate free surface position */
  //cal_fs(*topo.nx,*topo.ny,*topo.nz,stmp,*topo.fs);

  return;
}
