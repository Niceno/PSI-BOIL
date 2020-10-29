#include "commonheattransfer.h"

/***************************************************************************//*** 
*  \brief calculate ghost distance to interface resulting from heat trans resist
*******************************************************************************/
/* 
 * set of coords: cell with the interface, i.e. to
 * if cell_marker is > 0, material properties are of the other phase

table:

  .--------.-------.-------.
  | INT\CM |   +   |   -   |
  .--------.-------.-------.
  |   +    |   0   |   1   |
  .--------.-------.-------.
  |   -    |   1   |   0   |
  .--------.-------.-------.

  This convoluted way is used because I wanted to avoid checking if
  a cell, which is solid or in wall, is liquid/gas.
  Although update_at_walls should allow such a check, it is perhaps
  dangerous to rely on that. With the approach coded below, status of
  fluid cell is always checked.

*/

real CommonHeatTransfer::ghost_distance(const Comp & m, const Sign & cell_marker,
                                        const int i, const int j, const int k) 
                                        const {
  real res(0.);
  if(cell_marker*topo->sign_interface(i,j,k) < 0) {
    res = resistance_equivalent;
  }
  if(m==Comp::i()) {
    res /= std::max(boil::pico,fabs(topo->get_nx()[i][j][k]));
  } else if(m==Comp::j()) {
    res /= std::max(boil::pico,fabs(topo->get_ny()[i][j][k]));
  } else {
    res /= std::max(boil::pico,fabs(topo->get_nz()[i][j][k]));
  }
  return res;
}
