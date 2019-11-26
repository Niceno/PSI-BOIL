#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curvature() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  if       (bulk_curv_method==CurvMethod::HF()) {
    curv_HF();
  } else if(bulk_curv_method==CurvMethod::DivNorm()){
    if( !(phi.domain()->is_cartesian()) ) {
       boil::oout<<"Error: underdevelopment for curv_smooth for VOF on an "
                 <<" axisymmetric domain!"<<boil::endl;
       exit(0);
    }
    curv_smooth();
  } else if(bulk_curv_method==CurvMethod::none()) {
  } else {
    boil::oout<<"VOF::curvature: Curvature calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }

  return;
}

