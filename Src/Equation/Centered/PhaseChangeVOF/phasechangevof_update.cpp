#include "phasechangevof.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::update(const Scalar * diff_eddy) {

  boil::timer.start("phasechangevof update");

  /*------------+
  |  reset phi  |
  +------------*/
  initialize();

  /*----------------------+
  |  calculate mass flux  |
  +----------------------*/
  cal_massflux(diff_eddy);

  /*------------------------+
  |  calculate mass source  |
  +------------------------*/
  finalize();

  boil::timer.stop("phasechangevof update");

  return;
}

