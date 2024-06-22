#include "Include/psi-boil.h"

const real L  = 2.0;
const real D1 = 0.002;
const real D2 = 0.04;
const int  N  = 32;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /* stretched grid */
  Grid1D grid( Range<real>(-0.25*L, 0.75*L), 
               Range<real>(D1, D2), 
               N, 
               Periodic::no());

  //grid.plot("grid.eps");

  boil::timer.stop();
  boil::timer.report();
}
