#include "vof.h"

void VOF::curvature() {
  curvature(phi,true);

  return;
}

/******************************************************************************/
void VOF::curvature(const Scalar & Phi, const bool anc_flag) {
/***************************************************************************//**
*  \brief Calculate curvature
*******************************************************************************/

  if(curv_method==0) {
    curv_HF(Phi,anc_flag);
  } else {
    curv_smooth();
  }

  return;
}
