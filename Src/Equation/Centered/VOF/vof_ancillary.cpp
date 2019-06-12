#include "vof.h"

/******************************************************************************/
void VOF::ancillary() {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  /* with bnd update, the wall values of phi are distorted! 
   * update_at_walls should be called */
  phi.bnd_update();

#if 0
  /*---------------------------+
  |  flagging interface cells  |
  +---------------------------*/
  set_iflag();  /* done in curv_HF */
#endif

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm_cc(phi);

  boil::timer.start("vof ancillary");

  /* calculate alpha in cells */
  extract_alpha();

  /*  calculate free surface position */
#if 0
  cal_fs3();
#else
  cal_fs_interp();
#endif

  /* prerequisite for marching cubes */
  update_at_walls();

  /* for curvature calculations */
  norm_young(phi);
#if 0 /* norm_young produces the real-space normal vector as a byproduct */
  /* calculate the real-space normal vector */
  true_norm_vect(); 
#endif

  /* calculate area */
  cal_adens();
  //cal_adens_geom(adens);
  //set_adens(adensgeom);

  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();

  boil::timer.stop("vof ancillary");

  return;
}

