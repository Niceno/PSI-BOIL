#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::norm(const Scalar & color, const NormMethod & nm,
                     const bool extalp) {
/***************************************************************************//**
 \brief Calculate the normal vector and line constant, given color
    output: nx, ny, nz, nalpha
*******************************************************************************/

  if       (nm==NormMethod::ElviraXZ()) {
    /* elvira encapsulates nalpha calculations */
    norm_elvira(color,extalp);
  } else if(nm==NormMethod::Mixed()) {
    norm_mixed(color);
    if(extalp) extract_alpha(color);
  } else if(nm==NormMethod::Young()) {
    norm_young(color);
    if(extalp) extract_alpha(color);
  } else if(nm==NormMethod::CC()) {
    norm_cc(color);
    if(extalp) extract_alpha(color);
  } else {
    boil::oout<<"VOFaxisym::norm_axisymmetric: Normal vector calculation method "
              <<"not set properly! Exiting."<<boil::endl;
    exit(0);
  }

  //if(extalp) extract_alpha(color);

  return;
}


