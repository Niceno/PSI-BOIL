#include "vof.h"

/******************************************************************************/
void VOF::ancillary() {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  /* with bnd update, the wall values of phi are distorted! 
   * update_at_walls should be called */
  phi.bnd_update();

  boil::timer.start("vof ancillary");
  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
#if 1
  norm_mixed(phi);
  #if 0
  norm_cc(phi);
  #endif
#else
  norm_elvira(phi);
#endif

  /* iterate boundary normal vector */
  bdnorm();

  /* calculate alpha in cells */
  extract_alpha();

  /* prerequisite for marching cubes */
  update_at_walls();

  /* calculate the real-space normal vector */
  true_norm_vect(); 

  /*  calculate free surface position */
#if 1
  cal_fs3();
#else
  cal_fs_interp();
  if(!use_subgrid) fs_bnd_nosubgrid();
#endif

  /* calculate area */
  cal_adens();
#if 0
  cal_adens_geom(adens);
  set_adens(adensgeom);
#endif

  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();

  boil::timer.stop("vof ancillary");

  return;
}

