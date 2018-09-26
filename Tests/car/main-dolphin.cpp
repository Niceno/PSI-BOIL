#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  boil::plot = new PlotGMV();

  Grid1D gx( Range<real>(-150.00, 150.00), 128, Periodic::no() ); /* 6.0 */
  Grid1D gy( Range<real>(-450.00, 450.00), 256, Periodic::no() ); /* 3.0 */
  Grid1D gz( Range<real>(-160.00, 160.00), 128, Periodic::no() ); /* 1.5 */

//  Body car("car.stl");
  Body dolphin("dolphin.stl");

  Domain d(gx, gy, gz, dolphin);

//  boil::plot->plot(car, "car");
  boil::plot->plot(dolphin, "dolphin");

  boil::timer.stop();
  boil::timer.report();
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-dolphin.cpp,v 1.2 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
