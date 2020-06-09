#include "nucleation.h"

/******************************************************************************/
bool Nucleation::in_vapor(const int i, const int j, const int k) const {
/***************************************************************************//**
*  \brief test if in vapor taking into account sign
*******************************************************************************/
  if(sig==Sign::pos()) {
    return topo->below_interface((*clr)[i][j][k]);
  } else {
    return topo->above_interface((*clr)[i][j][k]);
  }
}

/******************************************************************************/
bool Nucleation::in_vapor(const real c) const {
/***************************************************************************//**
*  \brief test if in vapor taking into account sign
*******************************************************************************/
  if(sig==Sign::pos()) {
    return topo->below_interface(c);
  } else {
    return topo->above_interface(c);
  }
}
