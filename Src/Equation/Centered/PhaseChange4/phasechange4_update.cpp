#include "phasechange4.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange4::update(const Scalar * diff_eddy) {

  boil::timer.start("phasechange4 update");

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

  boil::timer.stop("phasechange4 update");

  return;
}

