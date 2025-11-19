#include "cavitypressure.h"

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
bool CavityPressure::interface(const Sign dir, const Comp m,
                               const int i, const int j, const int k) {

  return topo->interface(dir,m,i,j,k);
}
