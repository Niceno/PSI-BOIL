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
  //norm_cc(phi);
  norm_elvira(phi);

  /* calculate the real-space normal vector */
  true_norm_vect(); 

  boil::timer.start("vof ancillary");

  /* calculate alpha in cells */
  extract_alpha();

  /*  calculate free surface position */
#if 1
  cal_fs3();
#else
  cal_fs_interp();
#endif

  /* prerequisite for marching cubes */
  update_at_walls();

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

