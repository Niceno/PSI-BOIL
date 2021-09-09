#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::color_to_vf(Scalar & color, Scalar & vf, 
                           const bool nvec, const bool extalp, const bool bdn) {
/***************************************************************************//**
 \brief Solve the forward axisymmetric problem, i.e. calculate phi(alp,n,x),
    assuming the knowledge of Cartesian color
    output: vf = phi = volume fraction in an axisymmetric cell
            Ktmp = correction coefficient K = vf/c
*******************************************************************************/

  /* true = extract alpha */
  if(nvec) {
    norm(color,norm_method_advance,extalp);

    /* iterate boundary normal vector */
    if(bdn)
      bdnorm(color);
  }

  /* secondly, the forward problem is solved, Ktmp is overwritten */
  forward_axisymmetric(color,vf,Ktmp);

  return;
}


