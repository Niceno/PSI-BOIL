#include "tif.h"

/***************************************************************************//**
 *  Checks if the given cell is at an interface
******************************************************************************/
bool TIF::Interface(const int i, const int j, const int k) {
  if(Interface((*adens)[i][j][k])) {
    return true;
  } 
  return false;
}

bool TIF::Interface(const real heavi) {
  return heavi > boil::pico && heavi-1.0 < -boil::pico;
}
