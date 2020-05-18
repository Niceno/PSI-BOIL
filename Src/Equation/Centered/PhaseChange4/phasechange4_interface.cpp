#include "phasechange4.h"

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
bool PhaseChange4::interface(const Sign dir, const Comp m,
                             const int i, const int j, const int k) {

  return topo->interface(dir,m,i,j,k);
}

bool PhaseChange4::interface(const int i, const int j, const int k) {

  return topo->interface(i,j,k);
}
