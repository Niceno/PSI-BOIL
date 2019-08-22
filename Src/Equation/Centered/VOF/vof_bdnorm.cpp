#include "vof.h"

/******************************************************************************/
void VOF::bdnorm() {
/***************************************************************************//**
*  \brief iteration routine to obtain
          normal vector for cells adjacent to a wall or an immersed boundary
*         obtained using extrapolated volume fractions
*******************************************************************************/

  int niter = 2;
  for(int iter(0); iter<niter; ++iter) {

    /* calculate alpha in cells */
    extract_alpha();

    /* prerequisite for marching cubes */
    update_at_walls();
   
    /* calculate new normal vector near walls */
    extend_norm(phi);

  }

  return;
}
