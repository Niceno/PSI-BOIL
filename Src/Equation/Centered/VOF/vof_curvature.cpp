#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curvature() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  if(curv_method==0) {
    curv_HF();
  } else {
    curv_smooth();
  }

  return;
}

