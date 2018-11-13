#include "phasechangevof.h"

/***************************************************************************//**
 *  Returns interface temperature (calculated in EnthalpyTIF) 
*******************************************************************************/
real PhaseChangeVOF::Tint(const int i, const int j, const int k) {
  if (tif) return (*tif)[i][j][k];
  else return tsat;
}
     
