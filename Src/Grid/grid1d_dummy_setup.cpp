#include "grid1d.h"

void Grid1D::dummy_setup(const real dx) {

  allocate();

  /* distribute nodes inside */
  x_node[boil::BW]   = -dx/2.;
  x_node[boil::BW+1] = +dx/2.;

  correct_boundaries();

  return;
}
