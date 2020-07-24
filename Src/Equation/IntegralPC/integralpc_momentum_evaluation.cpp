#include "integralpc.h"

/******************************************************************************/
void IntegralPC::momentum_evaluation() {
/***************************************************************************//**
*  \brief solve momentum equation. 
*******************************************************************************/

  ns.discretize();
  ns.new_time_step(uvw_old);
  ns.grad(press);
  ns.solve(rr_ns);

  return;
}
