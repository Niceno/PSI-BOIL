#include "phasechange.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange::cal_massflux(const Scalar * diff_eddy) {

  /* calculate grad(tpr) */
  gradt(diff_eddy);

  /* calculate m */
  m(diff_eddy);
}

