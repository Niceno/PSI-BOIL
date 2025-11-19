#include "phasechangevof.h"

/***************************************************************************//**
 *  Returns interface temperature (calculated in Tifmodel) 
*******************************************************************************/
real PhaseChangeVOF::Tint(const int i, const int j, const int k) {
  return tifmodel.Tint(i,j,k);
}
     
