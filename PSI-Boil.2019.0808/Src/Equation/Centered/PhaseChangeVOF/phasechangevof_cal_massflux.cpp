#include "phasechangevof.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::cal_massflux(const Scalar * diff_eddy) {

  /* calculate grad(tpr) */
  gradt(diff_eddy);

  /* calculate m */
  m(diff_eddy);
}

