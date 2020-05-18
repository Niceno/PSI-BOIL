#include "enthalpyfd.h"

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
bool EnthalpyFD::interface(const Sign dir, const Comp m,
                           const int i, const int j, const int k) {

  return topo->interface(dir,m,i,j,k);
}

bool EnthalpyFD::interface_old(const Sign dir, const Comp m,
                               const int i, const int j, const int k) {

  return topo->interface_old(dir,m,i,j,k);
}
