#include "vof.h"

void VOF::ancillary() {
  ancillary(color());
  return;
}

/******************************************************************************/
void VOF::ancillary(Scalar & scp) {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  boil::timer.start("vof ancillary");

  /* reconstruct geometry */
  /* attention! PHI must be used here! */
  reconstruct_geometry(phi);

  if(topo_method==TopoMethod::Hybrid()) {

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

      heavi->cal_fs_interp(scp,fs,tol_wall,use_subgrid);

    }

    /* calculate area -> due to the Heaviside bindings, color is used, not scp */
    heavi->calculate_adens();

  } else if(topo_method==TopoMethod::Heaviside()) {

    /* under-developed part */
    bool use_subgrid;
    if       (subgrid_method == SubgridMethod::PLIC()) {
      use_subgrid = true;
    } else if(subgrid_method == SubgridMethod::None()) {
      use_subgrid = false;
    } else {
      boil::oout<<"VOF::ancillary: Underdevelopment!"
                <<" Exiting."<<boil::endl;
      exit(0);
    }

    heavi->topology(*topo,tol_wall,use_interp,use_subgrid);
  } else {
    boil::oout<<"VOF::ancillary: Topology calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }

  /* flag the interface */
  interfacial_flagging(scp);

  /* calculate scp in staggered cells */
  if(bndclr)
    cal_bndclr(scp);

  boil::timer.stop("vof ancillary");

  return;
}

