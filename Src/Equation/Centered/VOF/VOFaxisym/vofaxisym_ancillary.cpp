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

  /* reconstruct geometry */
  reconstruct_geometry(scp);

  /* calculate free surface position */
  if(!use_interp) {
    cal_fs3(scp);
  } else {
    cal_fs_interp(scp);
    if(!use_subgrid) 
      fs_bnd_nosubgrid(scp);
  }


  boil::timer.stop("vofaxisym ancillary");

  return;
}

