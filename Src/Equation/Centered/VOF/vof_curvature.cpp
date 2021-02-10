#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curvature() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  if (bulk_curv_method==CurvMethod::HF()) {
    curv_HF();
  } else if(bulk_curv_method==CurvMethod::none()) {
  } else {
    boil::oout<<"VOF::curvature: Curvature calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }

  return;
}

