#include "vof.h"

/******************************************************************************/
void VOF::bdnorm(Scalar & scp) {
/***************************************************************************//**
*  \brief iteration routine to obtain
          normal vector for cells adjacent to a wall or an immersed boundary
*         obtained using extrapolated volume fractions
*******************************************************************************/

  int niter = 2;
  for(int iter(0); iter<niter; ++iter) {

    /* prerequisite for marching cubes */
    if(update_at_walls_variable) {
      update_at_walls_custom(scp);
    } else {
      update_at_walls(scp);
    }
   
    /* calculate new normal vector near walls */
    normal_vector_near_bnd(scp,norm_method_advance);

    /* calculate alpha in cells */
    //extract_alpha(scp);
    extract_alpha_near_bnd(scp);
  }

  /* prerequisite for marching cubes */
  if(update_at_walls_variable) {
    update_at_walls_custom(scp);
  } else {
    update_at_walls(scp);
  }

  return;
}
