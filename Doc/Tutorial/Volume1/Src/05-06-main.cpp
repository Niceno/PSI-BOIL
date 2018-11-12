#include "Include/psi-boil.h"

const real LX = 4.0;  
const real LY = 1.0;
const int  NX =  64;
const int  NY =  16;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /* plot in Tecplot format */
  boil::plot = new PlotTEC();

  /* grids */
  Grid1D g_x( Range<real>(0,LX), NX, Periodic::no() );
  Grid1D g_y( Range<real>(0,LY), NY, Periodic::no() );

  /* create immeresed body */
  Body step("05-06-step.stl");

  /* plot step before immersion */
  boil::plot->plot(step, "step-before");

  /* create domain */
  Domain dom(g_x, g_y, g_y, &step);

  /* plot step after immersion */
  boil::plot->plot(step, "step-after");

  /* plot domain after immersion */
  boil::plot->plot(dom, "dom-after");

  boil::timer.stop();
  boil::timer.report();
}
