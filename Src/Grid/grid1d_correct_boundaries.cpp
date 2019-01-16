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

  for(int b=0; b<boil::BW; b++) {
    real dx1 = x_node[boil::BW + b + 1] - x_node[boil::BW];
    if(periodicN() == true) x_node[nc_in + boil::BW + b + 1] += dx1;

    real dxN = x_node[nc_in + boil::BW] - x_node[nc_in + boil::BW - b - 1];
    if(periodic1() == true) x_node[boil::BW - b - 1] -= dxN;
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
  if(periodic1()==true) dx_node[0]   = x_cell[nc_in+2*boil::BW] 
                                     - x_cell[nc_in+2*boil::BW-1];
  if(periodicN()==true) dx_node[nc_in+2*boil::BW] = x_cell[boil::BW+1] 
                                                  - x_cell[boil::BW];
}
