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

  /* with bnd update, the wall values of scp are distorted! 
   * update_at_walls should be called */
  scp.bnd_update();

  boil::timer.start("vof ancillary");
  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm(scp,norm_method_advance,true); /* alpha is extracted */

  /* iterate boundary normal vector */
  bdnorm(scp);

  /* prerequisite for marching cubes */
  update_at_walls(scp);

  /* calculate the real-space normal vector */
  true_norm_vect(nx,ny,nz,mx,my,mz); 

  /*  calculate free surface position */
#if 1
  cal_fs3(scp);
#else
  cal_fs_interp(scp);
  if(!use_subgrid) fs_bnd_nosubgrid(scp);
#endif

  /* calculate area -> due to the Heaviside bindings, phi is used, not scp */
  cal_adens();
#if 0
  cal_adens_geom(adens);
  set_adens(adensgeom);
#endif

  /* calculate scp in staggered cells */
  if(bndclr)
    cal_bndclr(scp);

  boil::timer.stop("vof ancillary");

  return;
}

