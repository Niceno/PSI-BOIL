#include "vofaxisym.h"

void VOFaxisym::ancillary() {
  ancillary(phi);
  return;
}

/******************************************************************************/
void VOFaxisym::ancillary(Scalar & scp) {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  boil::timer.start("vofaxisym ancillary");

  reconstruct_geometry(scp);

  boil::timer.stop("vofaxisym ancillary");

  return;
}

