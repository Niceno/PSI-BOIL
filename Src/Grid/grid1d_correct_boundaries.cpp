#include "grid1d.h"

/******************************************************************************/
void Grid1D::correct_boundaries() {

  /*------------------------+
  |                         |
  |  handle boundary nodes  |
  |                         |
  +------------------------*/
  for(int b=0; b<boil::BW; b++) {
    x_node[b]                        = x_node[boil::BW];
    x_node[nc_in + boil::BW + b + 1] = x_node[nc_in + boil::BW];
  }

  real dx1 = x_node[boil::BW + 1] - x_node[boil::BW];
  real dxN = x_node[nc_in + boil::BW] - x_node[nc_in + boil::BW - 1];

  for(int b=0; b<boil::BW; b++) {
    real dx1b = x_node[boil::BW + b + 1] - x_node[boil::BW];
    real dxNb = x_node[nc_in + boil::BW] - x_node[nc_in + boil::BW - b - 1];

    if(periodicN() == true) {
      x_node[nc_in + boil::BW + b + 1] += dx1b;
    } else {
      if     (cutoffN() == Cutoff::extrapolate())
        x_node[nc_in + boil::BW + b + 1] += real(b+1)*dxN;
      else if(cutoffN() == Cutoff::symmetry())
        x_node[nc_in + boil::BW + b + 1] += dxNb;
    }

    if(periodic1() == true) {
      x_node[boil::BW - b - 1] -= dxNb;
    } else {
      if     (cutoff1() == Cutoff::extrapolate())
        x_node[boil::BW - b - 1] -= real(b+1)*dx1;
      else if(cutoff1() == Cutoff::symmetry())
        x_node[boil::BW - b - 1] -= dx1b;
    }
  }

  /*--------------------------+            N = 4
  |                           |
  |  create cell coordinates  |    0   1   2   3   4   5
  |      and cell sizes       |  | o |-O-+-O-+-O-+-O-| o |
  |                           |  0   1   2   3   4   5   6
  +--------------------------*/
  for(int i=0; i<nc_in + 2*boil::BW; i++) {
     x_cell[i] = 0.5 * (x_node[i+1] + x_node[i]);
    dx_cell[i] =       (x_node[i+1] - x_node[i]);
  }

  /*--------------------------------------+            N = 4
  |                                       |
  |  create node (staggered cells) sizes  |    0   1   2   3   4   5
  |                                       |  | o |-O-+-O-+-O-+-O-| o |
  +---------------------------------------+  0   1   2   3   4   5   6  */
  for(int i=1; i<nc_in + 2*boil::BW; i++) 
    dx_node[i] = x_cell[i] - x_cell[i-1];

  dx_node[0]                  = 0;
  dx_node[nc_in + 2*boil::BW] = 0;
  if(periodic1()==true) {
    /* previously in error (indexing past the array + referring to wrong cell)
       dx_node[0] = x_cell[nc_in+2*boil::BW] - x_cell[nc_in+2*boil::BW-1]; */
    /* node[0] \equiv node[boil::BW+nc_in - boil::BW] = node[nc_in] */
    dx_node[0] = x_cell[nc_in] - x_cell[nc_in-1];
  } else {
    if     (cutoff1() == Cutoff::extrapolate())
      dx_node[0] = dx_node[1];
    else if(cutoff1() == Cutoff::symmetry())
      dx_node[0] = x_cell[2*boil::BW] - x_cell[2*boil::BW-1];
  }
  if(periodicN()==true) {
    /* previously in error (referring to wrong cell)
       dx_node[nc_in+2*boil::BW] = x_cell[boil::BW+1] - x_cell[boil::BW]; */
    dx_node[nc_in+2*boil::BW] = x_cell[2*boil::BW] - x_cell[2*boil::BW-1];
  } else {
    if     (cutoffN() == Cutoff::extrapolate())
      dx_node[nc_in+2*boil::BW] = dx_node[nc_in+2*boil::BW-1];
    else if(cutoffN() == Cutoff::symmetry())
      dx_node[nc_in+2*boil::BW] = x_cell[nc_in] - x_cell[nc_in-1];
  }

  return;
}
