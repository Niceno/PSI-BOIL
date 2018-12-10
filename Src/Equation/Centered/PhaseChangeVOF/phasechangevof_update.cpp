#include "phasechangevof.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::update(const Scalar * diff_eddy) {

  /*-----------------------------------------+
  |  calculate distance function + set flag  |
  +-----------------------------------------*/
  initialize();

  /*----------------------+
  |  calculate mass flux  |
  +----------------------*/
  cal_massflux(diff_eddy);

  /*------------------------+
  |  calculate mass source  |
  +------------------------*/
  finalize();
}
