#include "Include/psi-boil.h"

const real L = 2.0;
const int  N = 16;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /* non-periodic grid */
  Grid1D grid( Range<real>(-0.25*L, 0.75*L), N, Periodic::no());

  grid.print();
  grid.plot("grid.eps");

  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: 05-01-main.cpp,v 1.3 2008/11/17 19:23:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
