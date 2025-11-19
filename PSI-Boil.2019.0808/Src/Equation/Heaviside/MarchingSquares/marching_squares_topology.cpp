#include "marching_squares.h"

/******************************************************************************/
void MarchingSquares::topology(Topology & topo, const real tol_wall,
                               const bool use_interp, const bool use_subgrid) {
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

  (*topo.adens) = 0.0;
  (*nnx) = 0.0;
  (*nny) = 0.0;
  (*nzero) = 0.0;
  stmp = boil::unreal;

  for_vijk(stmp,i,j,k) {
    if(dom->ibody().off(i,j,k)) {
      continue;
    }

    /* lines get cleared with every call of this function */
    (*topo.adens)[i][j][k] = ad(i,j,k,lines);
 
#if 0
    /* more than single surface */
    if(extract_line_parameters(lines,
                               (*nnx)[i][j][k],(*nny)[i][j][k],
                               stmp[i][j][k]) > 1) {
      real nxval(0.0), nyval(0.0);
      if       (perpendicular==Comp::i()) {
        nxval = (nodalvals[i][j+1][k  ]+nodalvals[i][j+1][k+1])
               -(nodalvals[i][j  ][k  ]+nodalvals[i][j  ][k+1]);
        nyval = (nodalvals[i][j  ][k+1]+nodalvals[i][j+1][k+1])
               -(nodalvals[i][j  ][k  ]+nodalvals[i][j+1][k  ]);
      } else if(perpendicular==Comp::j()) {
        nxval = (nodalvals[i+1][j][k  ]+nodalvals[i+1][j][k+1])
               -(nodalvals[i  ][j][k  ]+nodalvals[i  ][j][k+1]);
        nyval = (nodalvals[i  ][j][k+1]+nodalvals[i+1][j][k+1])
               -(nodalvals[i  ][j][k  ]+nodalvals[i+1][j][k  ]);
      } else if(perpendicular==Comp::k()) {
        nxval = (nodalvals[i+1][j  ][k]+nodalvals[i+1][j+1][k])
               -(nodalvals[i  ][j  ][k]+nodalvals[i  ][j+1][k]);
        nyval = (nodalvals[i  ][j+1][k]+nodalvals[i+1][j+1][k])
               -(nodalvals[i  ][j  ][k]+nodalvals[i  ][j  ][k]);
      }

      real nsum = sqrt(nxval*nxval+nyval*nyval)+boil::pico;
      nxval /= nsum;
      nyval /= nsum;
    }
#else
    extract_line_parameters(lines,
                            (*nnx)[i][j][k],(*nny)[i][j][k],
                            stmp[i][j][k]);
#endif
  }

  nnx->exchange_all();
  nny->exchange_all();
  nzero->exchange_all();
  stmp.exchange_all();

  /* step two: calculate free surface position */
  if(!use_interp) {
    cal_fs_geom(*clr,*topo.nx,*topo.ny,*topo.nz,stmp,*topo.fs,tol_wall,use_subgrid);
  } else {
    cal_fs_interp(*clr,*topo.fs,tol_wall,use_subgrid);
  }

  return;
}
