#include "grid1d.h"

/******************************************************************************/
void Grid1D::correct_boundaries() {

  const real D1 = x_node[2]   - x_node[1];
  const real DN = x_node[N+1] - x_node[N];

  /*------------------------+
  |                         |
  |  handle boundary nodes  |
  |                         |
  +------------------------*/
  x_node[0]   = x_node[1];
  x_node[N+2] = x_node[N+1];
  if(periodic1()==true) x_node[0]   -= D1;
  if(periodicN()==true) x_node[N+2] += DN;

  /*--------------------------+            N = 4
  |                           |
  |  create cell coordinates  |    0   1   2   3   4   5
  |      and cell sizes       |  | o |-O-+-O-+-O-+-O-| o |
  |                           |  0   1   2   3   4   5   6    
  +--------------------------*/                                
  for(int i=0; i<N+2; i++) {
     x_cell[i] = 0.5 * (x_node[i+1] + x_node[i]);
    dx_cell[i] =       (x_node[i+1] - x_node[i]);
  }

  /*--------------------------------------+            N = 4
  |                                       |
  |  create node (staggered cells) sizes  |    0   1   2   3   4   5
  |                                       |  | o |-O-+-O-+-O-+-O-| o |
  +---------------------------------------+  0   1   2   3   4   5   6  */
  for(int i=1; i<N+2; i++) 
    dx_node[i] = x_cell[i] - x_cell[i-1];

  dx_node[0]   = 0;
  dx_node[N+2] = 0;
  if(periodic1()==true) dx_node[0]   = x_cell[N] - x_cell[N-1];
  if(periodicN()==true) dx_node[N+2] = x_cell[2] - x_cell[1];
}  
/*-----------------------------------------------------------------------------+
 '$Id: grid1d_correct_boundaries.cpp,v 1.3 2008/11/17 19:23:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
