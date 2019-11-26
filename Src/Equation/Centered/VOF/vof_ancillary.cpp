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

  if(topo_method==TopoMethod::Hybrid()) {

    /* calculate the real-space normal vector */
    true_norm_vect(nx,ny,nz,mx,my,mz); 

    /* calculate free surface position */
    if(!use_interp) {
      cal_fs3(scp);
    } else {
      heavi->cal_fs_interp(scp,fs);
      //cal_fs_interp(scp);
      if(!use_subgrid) 
        fs_bnd_nosubgrid(scp);
    }


    /* calculate area -> due to the Heaviside bindings, color is used, not scp */
    heavi->calculate_adens();

  } else if(topo_method==TopoMethod::Heaviside()) {
    heavi->topology(topo,use_interp);
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

