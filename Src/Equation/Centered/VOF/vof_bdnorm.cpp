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
    update_at_walls(scp);
   
    /* calculate new normal vector near walls */
    normal_vector_near_bnd(scp);

    /* calculate alpha in cells */
    //extract_alpha(scp);
    extract_alpha_near_bnd(scp);
  }

  return;
}
