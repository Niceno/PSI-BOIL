#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief Calculate grad(tpr) considering free surface position.
*******************************************************************************/

  /* calculate grad(tpr) */
  cal_gradt(diff_eddy);

#if 0
  boil::plot->plot(clr,txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  boil::plot->plot(clr,txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  exit(0);
#endif

#if 1
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
  /* extrapolate grad(tpr) */
  extrapolate(txv, 1);
  extrapolate(tyv, 1);
  extrapolate(tzv, 1);
  extrapolate(txl,-1);
  extrapolate(tyl,-1);
  extrapolate(tzl,-1);
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

#if 0
  /* subgrid treatment */
  subgrid();
#endif

#if 0
  boil::plot->plot(clr,tnl,tnv,"clr-tnl-tnv",time->current_step());
  exit(0);
#endif

  return;
}

