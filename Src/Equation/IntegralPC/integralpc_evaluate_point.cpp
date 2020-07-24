#include "integralpc.h"

/******************************************************************************/
real IntegralPC::evaluate_point(const real tval,const real pinf,
                                real & vpc_tpr, real & vpc_cav) {
/***************************************************************************//**
*  \brief evaluate discrepancy between the two phase change values.
*******************************************************************************/

  /* calculate cavity pressure value */
  real pcc = model(tval);

  /* restore temperature field */
  tpr = tpr_old;

  /* restore velocity field */
  for_m(m)
    uvw(m) = uvw_old(m);

  /* `standard' phase change */
  vpc_tpr = phase_change_evaluation(tval)*time.dt();

  /* cavity phase change */
  vpc_cav = cavity_evaluation(pinf,pcc);

  return vpc_cav-vpc_tpr;
}
