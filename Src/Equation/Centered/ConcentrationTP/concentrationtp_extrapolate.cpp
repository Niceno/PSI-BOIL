#include "concentrationtp.h"

/******************************************************************************/
void ConcentrationTP::extrapolate() {
/***************************************************************************//**
*  \brief Extrapolate eps across the interface.
*******************************************************************************/

  /* flagging (different from topology flag) */
  extrapolation_flag();

  /* extrapolate epsilon from vapor to liquid */
  topo->extrapolate(phi,-sig,{-sig},eflag);

  return;
}
