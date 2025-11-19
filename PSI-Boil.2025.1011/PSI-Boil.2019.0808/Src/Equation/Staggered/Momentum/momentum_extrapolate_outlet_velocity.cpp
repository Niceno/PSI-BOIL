#include "momentum.h"

/******************************************************************************/
void Momentum::extrapolate_outlet_velocity(const real fubo, const real ratio) {
  /*-------------------------------------------------+ 
  |  set range in which convective outlet will work  |
  +-------------------------------------------------*/
  Range<real> valid(0.97, 1.03);

  /*-------------------------+
  |  use convective outflow  |
  +-------------------------*/
  if( valid.contains(ratio) ) {
    boil::oout << "using convective outflow " << boil::endl;
    convective_outlet(u,fubo);

  /*-----------------------------------+
  |  use vanishing derivative outflow  |
  +-----------------------------------*/
  } else {
    boil::oout << "using vanishing derivative outflow " << boil::endl;
    vanishing_derivative_outlet(u);
  }

  return;
}
