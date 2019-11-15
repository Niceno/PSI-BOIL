#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::vf_to_color(const Scalar & vf, Scalar & color) {
/***************************************************************************//**
 \brief Solve the backward axisymmetric problem, i.e. calculate color(vf,n,x),
    assuming the knowledge of volumetric fraction and normal vector.
    information value.
    prerequisites: nx, ny, nz
    output: color = clr = area fraction in a Cartesian cell
            nalpha = line constant in a Cartesian cell
*******************************************************************************/

  /* firstly, alpha must be calculated */
  backward_axisymmetric(vf,nalpha);

  /* secondly, color is calculated with the forward Cartesian problem */
  //forward_cartesian(color);
  forward(color);

  return;
}


