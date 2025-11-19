#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::mass_flux(const Scalar * diff_eddy) {

  /* calculate heat flux */
  heat_flux(diff_eddy);

  /* calculate m */
  m();
}

