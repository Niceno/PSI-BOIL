#include "phasechangevof.h"

/******************************************************************************/
real PhaseChangeVOF::mdot_cut(real mdotval, real vfval, real & mcut) {
/***************************************************************************//**
*  \brief Cut off mdot to prevent over/undershoots of vf, preserving the cut val
*******************************************************************************/
  real dt=time->dt();
  real vfc = std::min(1.0,std::max(0.0,vfval));
  real mdotc;
  if(mdotval>=0.0) {
    mdotc = std::min(mdotval,rhol*vfc/dt);
    mcut = mdotval - mdotc;
  } else {
    mdotc = std::max(mdotval,-rhol*(1.0-vfc)/dt);
    mcut = mdotval - mdotc;
  }
  return mdotc;
}
