#include "phasechange4.h"

/******************************************************************************/
real PhaseChange4::mdot_cut(real mdotval, real vfval, real & mcut) {
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

/******************************************************************************/
real PhaseChange4::vfs_cut(real vfsval, real vfval) {
/***************************************************************************//**
*  \brief Cut off vfs to prevent over/undershoots of vf
*******************************************************************************/
  real dt=time->dt();
  real vfc = std::min(1.0,std::max(0.0,vfval));
  real vfsc;
  if(vfsval<=0.0) {
    vfsc = -std::min(-vfsval,vfc/dt);
  } else {
    vfsc = std::min(vfsval,(1.0-vfc)/dt);
  }
  return vfsc;
}
