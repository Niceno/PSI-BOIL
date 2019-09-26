#include "vofaxisym.h"

/******************************************************************************/
real VOFaxisym::linf_scalar_error(const Scalar & sca, const Scalar & scb) {
/***************************************************************************//**
 \brief Calculate the Linf difference between two scalar fields.
    output: Linf
*******************************************************************************/

  real linferr(-1.0);
  for_vijk(sca,i,j,k) {
    if(dom->ibody().off(i,j,k)) continue;
    real sdiff = fabs(sca[i][j][k] - scb[i][j][k]);
    if(sdiff>linferr)
      linferr = sdiff;
  }
  boil::cart.max_real(&linferr);

  return linferr;
}
