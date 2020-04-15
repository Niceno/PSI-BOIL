#include "phasechange4.h"

/***************************************************************************//**
 *  Returns interface temperature (calculated in Tifmodel) 
*******************************************************************************/
real PhaseChange4::Tint(const int i, const int j, const int k) {
  return tifmodel.Tint(i,j,k);
}
     
