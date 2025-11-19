#include "vof.h"

/******************************************************************************/
void VOF::norm(const Scalar & color, const NormMethod & nm,
               const bool extalp) {
/***************************************************************************//**
 \brief Calculate the normal vector and line constant, given color
    output: nx, ny, nz, nalpha
*******************************************************************************/

  if       (nm==NormMethod::Mixed()) {
    norm_mixed(color);
    if(extalp) extract_alpha(color);
  } else if(nm==NormMethod::Young()) {
    norm_young(color);
    if(extalp) extract_alpha(color);
  } else if(nm==NormMethod::CC()) {
    norm_cc(color);
    if(extalp) extract_alpha(color);
  } else if(nm==NormMethod::ElviraXZ()) {
    norm_elvira(color,extalp);
  } else if(nm==NormMethod::ElviraXY()) {
    norm_elvira(color,extalp);
  } else if(nm==NormMethod::ElviraYZ()) {
    norm_elvira(color,extalp);
  } else {
    /* default */
    boil::oout<<"VOF::norm: Normal vector calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }

  //if(extalp) extract_alpha(color);

  return;
}


