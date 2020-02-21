#include "concentrationtp.h"

/******************************************************************************/
void ConcentrationTP::extrapolate() {
/***************************************************************************//**
*  \brief Extrapolate eps across the interface.
*******************************************************************************/

  /* flagging (different from topology flag) */
  extrapolation_flag();

  /* extrapolate epsilon from vapor to liquid */
  topo->extrapolate(phi,Sign::pos(),eflag);

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: phasechange_normalize.cpp,v 1.1 2014/02/03 14:28:07 sato Exp $'/
+-----------------------------------------------------------------------------*/
