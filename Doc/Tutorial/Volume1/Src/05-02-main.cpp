#include "Include/psi-boil.h"

const real L = 2.0;
const int  N = 16;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /* periodic grid */
  Grid1D grid( Range<real>(-0.25*L, 0.75*L), N, Periodic::yes());

  grid.print();
  grid.plot("grid.eps");

  boil::timer.stop();
  boil::timer.report();
}
