#include "vof.h"

/******************************************************************************/
void VOF::ancillary() {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  /*--------------------------+
  |  flagging interface cells  |
  +--------------------------*/
  set_iflag();

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm_cc(phi);

  /* calculate alpha in cells */
  extract_alpha();

  /*  calculate free surface position */
  cal_fs3();

  /* prerequisite for marching cubes */
  update_at_walls();

  /* calculate area */
  cal_adens();

  /* calculate the real-space normal vector */
  true_norm_vect(); 

  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();

  return;
}

