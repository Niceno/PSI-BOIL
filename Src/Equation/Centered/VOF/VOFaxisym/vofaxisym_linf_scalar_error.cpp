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
    real sdiff = fabs(std::max(0.0,std::min(1.0,sca[i][j][k]))
                     -std::max(0.0,std::min(1.0,scb[i][j][k])));
    if(sdiff>linferr)
      linferr = sdiff;
  }
  boil::cart.max_real(&linferr);

  return linferr;
}
