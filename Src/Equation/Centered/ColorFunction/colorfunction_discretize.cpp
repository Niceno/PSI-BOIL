#include "colorfunction.h"

/******************************************************************************/
void ColorFunction::discretize() {

  boil::timer.start("colorfunction discretize");

  if(phi.dxc(1) != phi.dyc(1) ||
     phi.dxc(1) != phi.dzc(1) ||  
     phi.dyc(1) != phi.dzc(1)) 
    boil::oout << "Warning: the mesh is not uniform in all directions!" 
               << boil::endl 
               << "         It is not good for surface tracking with LS."
               << boil::endl;

  /*----------------------------------------------------------+
  |  compute local conductivity                               |
  |                                                           |
  |  smaller lcnd - narrower band, higher parasitic currents  |
  +----------------------------------------------------------*/
  real delta = boil::minr(phi.dxc(1), phi.dyc(1), phi.dzc(1));
  OPR(delta);

  const real lambda = 0.5 * delta;

  v_fluid.lambda( lambda );

  ldt = 0.15 * delta*delta / lambda;  OPR(ldt);

  /* create only diffusion matrix ... */
  create_system_diffusive(flu->lambda());

  /* ... and correct on the boundaries */
  create_system_bnd();   

  boil::timer.stop("colorfunction discretize");
}
