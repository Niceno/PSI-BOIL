#include "phasechange.h"
using namespace std;

/******************************************************************************/
void PhaseChange::gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief Calculate grad(tpr) considering free surface position.
*******************************************************************************/

  /* calculate grad(tpr) */
  cal_gradt(diff_eddy);
  
  /* set iflag */
  setflag();

  /* extrapolate grad(tpr) */
  ext_gradt(txv, 1);
  ext_gradt(tyv, 1);
  ext_gradt(tzv, 1);
  ext_gradt(txl,-1);
  ext_gradt(tyl,-1);
  ext_gradt(tzl,-1);

  //boil::plot->plot(clr,txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  //boil::plot->plot(clr,txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  //exit(0);

  return;
}
