#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::mdot_cut() {
/***************************************************************************//**
*  \brief Cut off mdot
*******************************************************************************/
  real dt=time->dt();
  for_vijk(tpr,i,j,k){
    real mdotc=phi[i][j][k];
    real clrc =clr[i][j][k];
    clrc = std::min(1.0,std::max(0.0,clrc));
    mdotc = std::min(mdotc,rhol*clrc/dt);
    phi[i][j][k] = std::max(mdotc,-rhol*(1.0-clrc)/dt);
  }
  phi.exchange_all();

  return;
}
/******************************************************************************/
real PhaseChangeVOF::mdot_cut(real mdotc, real clrc) {
/***************************************************************************//**
*  \brief Cut off mdot
*******************************************************************************/
  real dt=time->dt();
  clrc = std::min(1.0,std::max(0.0,clrc));
  mdotc = std::min(mdotc,rhol*clrc/dt);
  return std::max(mdotc,-rhol*(1.0-clrc)/dt);
}
