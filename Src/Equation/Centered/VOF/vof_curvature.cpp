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
    if( !(phi.domain()->is_cartesian()) ) {
       boil::oout<<"Error: underdevelopment for curv_smooth for VOF on an "
                 <<" axisymmetric domain!"<<boil::endl;
       exit(0);
    }
    curv_smooth();
  }

  return;
}

