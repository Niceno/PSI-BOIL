#include "cavitypressure.h"

/***************************************************************************//**
*  Corrects for the inactive gas cells. 
*******************************************************************************/
real CavityPressure::update_rhs() {

  return Pressure::update_rhs_pressure();
}
