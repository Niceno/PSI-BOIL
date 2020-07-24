#include "integralpc.h"

/******************************************************************************/
real IntegralPC::cavity_evaluation(const real pinf, const real pcc) {
/***************************************************************************//**
*  \brief evaluate volumetric phase change due to over-pressure.
*******************************************************************************/

  /* step one: reset outflow velocity */
  for_m(m)
    for_avmijk(uvw_cav,m,i,j,k)
      uvw_cav[m][i][j][k]=0.0;

  /* step two: set pressure */
  capr.set_cavity_pressure(pcc+pinf);

  /* step three: solve cavity pressure */
  press = pinf;
  capr.discretize();
  capr.coarsen();
  if(multigrid_cavity.cycle(c0,c1,rr_cav,mi))
    OMS(converged);
  press.exchange();

  /* step four: solve momentum equation */
  momentum_evaluation();

  /* step five: correct boundaries for velocity */
  ns.vanishing_derivative_outlet(uvw);

  /* step six: calculate outlet */
  return -uvw.bnd_flow(BndType::outlet()); /* m3/s */
}
