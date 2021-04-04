#include "nucleation.h"

/******************************************************************************/
bool Nucleation::in_vapor(const int i, const int j, const int k) const {
/***************************************************************************//**
*  \brief test if in vapor taking into account sign
*******************************************************************************/
  if(matter_sig==Sign::pos()) {
    return cht->topo->below_interface(i,j,k);
  } else {
    return cht->topo->above_interface(i,j,k);
  }
}

/******************************************************************************/
bool Nucleation::in_vapor(const real c) const {
/***************************************************************************//**
*  \brief test if in vapor taking into account sign
*******************************************************************************/
  if(matter_sig==Sign::pos()) {
    return cht->topo->below_interface(c);
  } else {
    return cht->topo->above_interface(c);
  }
}
