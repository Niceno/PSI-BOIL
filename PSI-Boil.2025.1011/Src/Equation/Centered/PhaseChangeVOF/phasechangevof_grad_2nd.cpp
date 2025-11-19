#include "phasechangevof.h"

/******************************************************************************/
real PhaseChangeVOF::grad_2nd(const real tm0, const real tm1, const real tm2,
                              const real dm1, const real dm2){
/***************************************************************************//**
*  \brief calculate second order accurate grad 
*******************************************************************************/
  real denom = dm1*dm2*(dm2-dm1);

  return dm2*dm2*(tm1-tm0)/denom - dm1*dm1*(tm2-tm0)/denom;
}
                
