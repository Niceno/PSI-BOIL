#include "Include/psi-boil.h"

const real L = 3.2;
const real D = 0.01;
const int  N = 32;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /* uniform grid */
  Grid1D g_uni( Range<real>(-0.5*L, 0.5*L), N, Periodic::no() );

  /* stretched grid */
  Grid1D g_str( Range<real>(0,L), Range<real>(D,D), N*2, Periodic::no());

  /* create domain */
  Domain dom(g_uni, g_uni, g_str);

  /* plot the domain */
  boil::plot = new PlotTEC();
  boil::plot->plot(dom, "dom");

  boil::timer.stop();
  boil::timer.report();
}
