#include "marching_squares.h"

/******************************************************************************/
void MarchingSquares::evaluate_nodes() {
/***************************************************************************//**
*  \brief Evaluate values of color at node points. Staggering is decided based
*         on offsets introduced by the perpendicular direction
*         Results: nodalvals
*******************************************************************************/

  clr->exchange_all();

  for_vijk(nodalvals,i,j,k) {
    real nn(0.0);
    for(int idx=ofx;idx<1;idx++)
      for(int jdx=ofy;jdx<1;jdx++)
        for(int kdx=ofz;kdx<1;kdx++)
          nn+=std::max(0.0,std::min(1.0,(*clr)[i+idx][j+jdx][k+kdx]));
    nn /= 4.0;
    nodalvals[i][j][k] = nn;
  }

  //nodalvals.exchange_all();

  return;
}
