#include "vof.h"

void VOF::reconstruct_geometry() {
  reconstruct_geometry(phi);
  return;
}

/******************************************************************************/
void VOF::reconstruct_geometry(Scalar & scp) {
/***************************************************************************//**
 \brief Reconstruct the geometry of the interface, i.e. calculate the normal
        vector, line constant alpha.
    line: vm1*x + vm2*z = alpha
    output: nx, ny = 0.0, nz, nalpha, c at walls
*******************************************************************************/

  /* normal vector at cell center */
  norm(scp,norm_method_advance,true); /* alpha is extracted */

  /* iterate boundary normal vector */
  bdnorm(scp);

  /* prerequisite for marching cubes */
  update_at_walls(scp);

  return;
}
