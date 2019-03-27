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
  ext_gradt(tnv, 1);
  ext_gradt(tnl,-1);

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
    for_avi(tzv,i)
       boil::oout<<i<<" "<<clr[i][48][48]<<" "<<txl[i][48][48]<<" "<<txv[i][48][48]
                 <<" | "<<clr[i][48][49]<<" "<<txl[i][48][49]<<" "<<txv[i][48][49]
                 <<" | "<<clr[i][49][48]<<" "<<txl[i][49][48]<<" "<<txv[i][49][48]
                 <<" | "<<clr[i][49][49]<<" "<<txl[i][49][49]<<" "<<txv[i][49][49]<<boil::endl;
#endif

#if 0
  //boil::plot->plot(clr,txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  boil::plot->plot(clr,tnl,tnv,"clr-tnl-tnv",time->current_step());
  exit(0);
#endif

  return;
}

