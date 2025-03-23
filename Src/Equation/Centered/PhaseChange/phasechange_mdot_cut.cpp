#include "phasechange.h"
using namespace boil;

/******************************************************************************/
void PhaseChange::mdot_cut() {
/***************************************************************************//**
*  \brief Cut off mdot
*******************************************************************************/
  real dt=time->dt();
  for_vijk(tpr,i,j,k){
    real mdotc=phi[i][j][k];
    real clrc =clr[i][j][k];
    clrc = minr(1.0,maxr(0.0,clrc));
    mdotc = minr(mdotc,rhol*clrc/dt);
    phi[i][j][k] = maxr(mdotc,-rhol*(1.0-clrc)/dt);
  }
  phi.exchange_all();

  return;
}
/******************************************************************************/
real PhaseChange::mdot_cut(real mdotc, real clrc) {
/***************************************************************************//**
*  \brief Cut off mdot
*******************************************************************************/
  real dt=time->dt();
  clrc = minr(1.0,maxr(0.0,clrc));
  mdotc = minr(mdotc,rhol*clrc/dt);
  return maxr(mdotc,-rhol*(1.0-clrc)/dt);
}
