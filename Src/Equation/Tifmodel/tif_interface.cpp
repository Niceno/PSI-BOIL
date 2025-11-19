#include "tif.h"

/***************************************************************************//**
 *  Checks if the given cell is at an interface
******************************************************************************/
bool TIF::Interface(const int i, const int j, const int k) {
  if(adens[i][j][k]>boil::pico) {
    return true;
  } 

  return false;
}
