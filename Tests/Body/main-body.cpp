#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  boil::plot = new PlotGMV();

  Grid1D gx( Range<real>(-1.5001, 1.5001),  50, Periodic::no() ); 
  Grid1D gy( Range<real>(-0.7501, 0.7501),  25, Periodic::no() ); 

//bad:  Grid1D gx( Range<real>(-1.50, 1.50),  60, Periodic::no() ); 
//bad:  Grid1D gy( Range<real>(-0.75, 0.75),  30, Periodic::no() ); 

  Body complet("sphere_and_cube.stl");
//  Body complet("psi.stl");

  Domain d(gx, gy, gy, & complet);

  boil::plot->plot(complet, "sphere_and_cube");

  boil::timer.stop();
  boil::timer.report();
}	
