#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::subgrid() {
/***************************************************************************//**
*  \brief resolve the subgrid model of near-wall heat transfer 
*******************************************************************************/
  
  /* flag subgrid cells */
  subgrid_setflag();

  //boil::plot->plot(clr,txl,tyl,tzl, "ib+I+clr-tx-ty-tz", time->current_step());

#if 0 /* currently, a simpler method is implemented, assuming 1D gradient */
  /* extrapolate tx,ty,tz to flagged cells */
  extrapolate(txv, 1);
  extrapolate(tyv, 1);
  extrapolate(tzv, 1);
  extrapolate(txl,-1);
  extrapolate(tyl,-1);
  extrapolate(tzl,-1);
#endif

  /* correct tnl and tnv in subgrid cells */
  subgrid_gradt();

  return;
}
