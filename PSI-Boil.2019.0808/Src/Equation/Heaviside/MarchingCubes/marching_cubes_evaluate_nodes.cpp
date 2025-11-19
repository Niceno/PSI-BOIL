#include "marching_cubes.h"

/******************************************************************************/
void MarchingCubes::evaluate_nodes() {
/***************************************************************************//**
*  \brief Evaluate values of color at node points.
*         Results: nodalvals
*******************************************************************************/

  clr->exchange_all();

  for_vijk(nodalvals,i,j,k) {
    real nn(0.0);
    for(int idx=-1;idx<1;idx++)
      for(int jdx=-1;jdx<1;jdx++)
        for(int kdx=-1;kdx<1;kdx++)
          nn+=std::max(0.0,std::min(1.0,(*clr)[i+idx][j+jdx][k+kdx]));
    nn /= 8.0;
    nodalvals[i][j][k] = nn;
  }

  //nodalvals.exchange_all();

  return;
}
