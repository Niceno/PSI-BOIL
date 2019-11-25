#include "vof.h"

void VOF::ancillary() {
  ancillary(phi);
  return;
}

/******************************************************************************/
void VOF::ancillary(Scalar & scp) {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  boil::timer.start("vof ancillary");

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

  /* calculate area -> due to the Heaviside bindings, color is used, not scp */
  heavi->calculate_adens();

  /* flag the interface */
  interfacial_flagging(scp);

  /* calculate scp in staggered cells */
  if(bndclr)
    cal_bndclr(scp);

  boil::timer.stop("vof ancillary");

  return;
}

