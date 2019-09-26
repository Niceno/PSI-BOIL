#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::norm_axisymmetric(const Scalar & color) {
/***************************************************************************//**
 \brief Calculate the normal vector and line constant, given color
    output: nx, ny, nz, nalpha
*******************************************************************************/

  /* firstly, alpha and nx,ny,nz must be calculated */
  if        (norm_method_advance==NormMethod::ElviraXZ()) {
    /* elvira encapsulates nalpha calculations */
    norm_elvira(color);
  } else if(norm_method_advance==NormMethod::Mixed()) {
    norm_mixed(color);
    extract_alpha(color);
  } else if(norm_method_advance==NormMethod::Young()) {
    norm_young(color);
    extract_alpha(color);
  } else if(norm_method_advance==NormMethod::CC()) {
    norm_cc(color);
    extract_alpha(color);
  } else {
    boil::aout<<"VOFaxisym::reconstruct_geom: Normal vector calculation method "
              <<"not set properly! Exiting."<<boil::endl;
    exit(0);
  }

  //extract_alpha(color);

  return;
}


