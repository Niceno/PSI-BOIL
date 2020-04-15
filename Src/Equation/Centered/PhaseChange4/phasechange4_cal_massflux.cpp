#include "phasechange4.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange4::cal_massflux(const Scalar * diff_eddy) {

  /* calculate grad(tpr) */
  gradt(diff_eddy);

  /* calculate m */
  m(diff_eddy);
}

