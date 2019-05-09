#include "vof.h"

/******************************************************************************/
void VOF::ancillary() {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

#if 1
  /*--------------------------+
  |  flagging interface cells  |
  +--------------------------*/
  set_iflag();  /* done in curv_HF */

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  //norm_cc(phi);
  norm_young(phi);
#else
  curv_HF();
#endif
  boil::timer.start("vof ancillary");

  /* calculate alpha in cells */
  extract_alpha();

  /*  calculate free surface position */
  //cal_fs3();
  cal_fs_interp();

  /* prerequisite for marching cubes */
  update_at_walls();

  /* calculate the real-space normal vector */
  true_norm_vect(); 

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

