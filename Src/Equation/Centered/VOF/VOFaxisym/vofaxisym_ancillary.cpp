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


  /* calculate the real-space normal vector */
  true_norm_vect(nx,ny,nz,mx,my,mz); 

  /* calculate free surface position */
  if(!use_interp) {
    cal_fs3(scp);
  } else {
    cal_fs_interp(scp);
    if(!use_subgrid) 
      fs_bnd_nosubgrid(scp);
  }

  /* flag the interface */
  interfacial_flagging(scp);

  boil::timer.stop("vofaxisym ancillary");

  return;
}

