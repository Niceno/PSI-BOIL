#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief Calculate grad(tpr) considering free surface position.
*******************************************************************************/

  /* calculate grad(tpr) */
  cal_gradt(diff_eddy);

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
  ext_gradt(txv, 1, Comp::i());
  ext_gradt(tyv, 1, Comp::j());
  ext_gradt(tzv, 1, Comp::k());
  ext_gradt(txl,-1, Comp::i());
  ext_gradt(tyl,-1, Comp::j());
  ext_gradt(tzl,-1, Comp::k());
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

