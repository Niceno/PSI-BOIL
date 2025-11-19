#include "vof.h"

void VOF::ancillary(const bool reconstruct) {
  ancillary(color(),reconstruct);
  return;
}

/******************************************************************************/
void VOF::ancillary(Scalar & scp, const bool reconstruct) {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  boil::timer.start("vof ancillary");

  /* reconstruct geometry */
  /* attention! PHI must be used here! */
  if(reconstruct)
    reconstruct_geometry(phi);

  /* calculate the real-space normal vector */
  true_norm_vect(nx,ny,nz,mx,my,mz); 

  /* calculate free surface position */
  if(!use_interp) {
    cal_fs3(scp);
  } else {

    /* with interpolation, only binary subgrid division is possible */
    bool use_subgrid;
    if(subgrid_method == SubgridMethod::None()) {
      use_subgrid = false;
    } else {
      use_subgrid = true;
    }

    topo->cal_fs_interp(scp,fs,tol_wall,use_subgrid);

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

