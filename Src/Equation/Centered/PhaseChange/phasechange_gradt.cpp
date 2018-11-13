#include "phasechange.h"
using namespace std;

/******************************************************************************/
void PhaseChange::gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief Calculate grad(tpr) considering free surface position.
*******************************************************************************/

  /* calculate grad(tpr) */
  cal_gradt(diff_eddy);

  /* extrapolate grad(tpr) */
  ext_gradt(txv, 1);
  ext_gradt(tyv, 1);
  ext_gradt(tzv, 1);
  ext_gradt(txl,-1);
  ext_gradt(tyl,-1);
  ext_gradt(tzl,-1);

#if 0 
   if(time->current_step()==200) {
     for_avik(tzv,i,k)
       boil::oout<<i<<" "<<k<<" "<<clr[i][1][k]<<" "<<tzl[i][1][k]<<" "<<tzv[i][1][k]<<boil::endl;
     exit(0); 
   }
#endif

  //boil::plot->plot(clr,txv,tyv,tzv,"clr-txv-tyv-tzv",time->current_step());
  //boil::plot->plot(clr,txl,tyl,tzl,"clr-txl-tyl-tzl",time->current_step());
  //exit(0);

  return;
}

