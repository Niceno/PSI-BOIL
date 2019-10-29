#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::color_to_vf(const Scalar & color, Scalar & vf, const bool extalp) {
/***************************************************************************//**
 \brief Solve the forward axisymmetric problem, i.e. calculate phi(alp,n,x),
    assuming the knowledge of Cartesian color
    output: vf = phi = volume fraction in an axisymmetric cell
            Ktmp = correction coefficient K = vf/c
*******************************************************************************/

  /* true = extract alpha */
  norm(color,norm_method_advance,extalp);

  /* secondly, the forward problem is solved, Ktmp is overwritten */
  forward_axisymmetric(color,vf,Ktmp);

  return;
}


