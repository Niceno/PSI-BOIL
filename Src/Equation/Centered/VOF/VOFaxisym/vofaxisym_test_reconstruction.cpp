#include "vofaxisym.h"

/******************************************************************************/
real VOFaxisym::test_reconstruction(const Scalar & color, const Scalar & vf) {
/***************************************************************************//**
 \brief Calculate the error of reconstruction by solving the forward problem
*******************************************************************************/

  /* Ktmp, nx, ny, nz, nalpha are internally overwritten in this function! */
  color_to_vf(color,axistmp);

  /* calculate Linf error norm of reconstruction in phi-space */
  real linferr = linf_scalar_error(axistmp,vf);

  boil::oout<<"VOFaxisym::test_reconstruction: error= "<<time->current_time()
           <<" "<<linferr<<boil::endl;

  return linferr;
}