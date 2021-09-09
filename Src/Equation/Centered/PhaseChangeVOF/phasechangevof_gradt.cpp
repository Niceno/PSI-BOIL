#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief Calculate grad(tpr) considering free surface position.
*******************************************************************************/

  /* used by gradt_ib and gradt8 */
  /* -> increased precision near cell centre */
  /* achieved by switching to upwind/downwind differencing */
  /* it is also possible to use purely upw/dwnw through the upwind_flag */
  /* this is however not recommended <- lower precision */
  calculate_node_temperature(diff_eddy);

  /* calculate grad(tpr) in fluid */
  cal_gradt_fluid(diff_eddy);

  /* correct grad(tpr) in fluid near ib */
  correct_gradt_at_ib(diff_eddy);

  /* calculate grad(tpr) in ib */
  cal_gradt_ib(diff_eddy);
  
  txl.bnd_update();
  tyl.bnd_update();
  tzl.bnd_update();
  txv.bnd_update();
  tyv.bnd_update();
  tzv.bnd_update();

  /* insert boundary condition at walls */
  insert_bc_gradt_at_walls(diff_eddy);

  txl.exchange_all();
  tyl.exchange_all();
  tzl.exchange_all();
  txv.exchange_all();
  tyv.exchange_all();
  tzv.exchange_all();

#if 0
  boil::plot->plot(clr,txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  boil::plot->plot(clr,txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  exit(0);
#endif

#if 0 /* since values outside the fluid domain are used, extrapolation must 
         come before the scalar product */
  /* calculate the normal component */
  for_aijk(i,j,k) {
    tnv[i][j][k] = txv[i][j][k]*nx[i][j][k]
                 + tyv[i][j][k]*ny[i][j][k]
                 + tzv[i][j][k]*nz[i][j][k];
    tnl[i][j][k] = txl[i][j][k]*nx[i][j][k]
                 + tyl[i][j][k]*ny[i][j][k]
                 + tzl[i][j][k]*nz[i][j][k];
  }
  /* extrapolate grad(tpr) */
  extrapolate(tnv, 1);
  extrapolate(tnl,-1);
#else
  /* extrapolate grad(tpr), values in brackets indicate
     iflag values for extrapolated cells */
  topo->extrapolate(txv,Sign::pos(),{1,2});
  topo->extrapolate(tyv,Sign::pos(),{1,2});
  topo->extrapolate(tzv,Sign::pos(),{1,2});
  topo->extrapolate(txl,Sign::neg(),{-1,-2});
  topo->extrapolate(tyl,Sign::neg(),{-1,-2});
  topo->extrapolate(tzl,Sign::neg(),{-1,-2});
  /* calculate the normal component */
  for_aijk(i,j,k) {
    tnv[i][j][k] = txv[i][j][k]*nx[i][j][k]
                 + tyv[i][j][k]*ny[i][j][k]
                 + tzv[i][j][k]*nz[i][j][k];
    tnl[i][j][k] = txl[i][j][k]*nx[i][j][k]
                 + tyl[i][j][k]*ny[i][j][k]
                 + tzl[i][j][k]*nz[i][j][k];
  }
  
#if 0
  boil::plot->plot(clr,txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  boil::plot->plot(clr,txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  exit(0);
#endif

#endif
  
  /* for increased stability in the near-wall region */
  if(near_wall_modelling)
    near_wall_model(diff_eddy);

#if 0
  boil::plot->plot(clr,tnl,tnv,"clr-tnl-tnv",time->current_step());
  exit(0);
#endif

  return;
}

