#include "integralpc.h"

/******************************************************************************/
real IntegralPC::phase_change_evaluation(const real tval) {
/***************************************************************************//**
*  \brief evaluate volumetric phase change due to heat transfer.
*******************************************************************************/

  /* set interfacial temperature field */
  tsat.set_tref(tval);

  /* solve enthalpy equation */
  enthFD.discretize();
  enthFD.new_time_step();
  enthFD.solve(ResRat(1e-16),"enthFD");

  /* solve phase change */
  f = 0.0;
  pc.update();

  /* sum up */
  real vpc_tpr(0.0);
  for_vijk(f,i,j,k)
    vpc_tpr+=f[i][j][k];

  boil::cart.sum_real(&vpc_tpr);

  return vpc_tpr;
}
