#include "phasechange.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange::initialize() {

  phi=0.0;
  /*------------------------------+
  |  calculate distance function  |
  +------------------------------*/
  /* calculate distance function */
  distfunc(clr,24);

  /* calculate normal vectors from distance function */
  gradphic(dist);
#ifdef DEBUG
  boil::plot->plot(clr, tpr, dist, "clr-tpr-dist",  time->current_step());
  boil::plot->plot(nx, ny, nz, "nx-ny-nz",  time->current_step());
#endif
  
  /* set iflag */
  setflag();
}

