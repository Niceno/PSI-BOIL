#include "cavitypressure.h"

/******************************************************************************/
CavityPressure::CavityPressure(const Scalar & PHI,
                               const Scalar & F,
                               const Vector & U,
                               Times & T,
                               Linear * sm,
                               Matter * liq,
                               Topology * TOPO,
                               const Property * TENS,
                               const Scalar & CURV,
                               Sign SIG) :
  /* call parent's constructor. NULL is for solid */
  Centered( PHI.domain(), PHI, F, & U, T, liq, NULL, sm ),
  topo(TOPO),
  fs(TOPO->fs),
  iflag(TOPO->iflag),
  matter_sig(SIG),
  sigma(TENS),
  kappa(&CURV)
{
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == sm->domain());

  cavity_pressure = 0.0;

  /* the in_gas function is initialized using a lambda expression */
  /* clr = 1 is liquid */
  if(SIG>0) {
    in_gas = [this](const int i, const int j, const int k)
                   { return topo->below_interface(i,j,k); };
  /* clr = 0 is liquid */
  } else {
    in_gas = [this](const int i, const int j, const int k)
                   { return topo->above_interface(i,j,k); };
  }

  phi.bnd_update();
  discretize();
}	

/******************************************************************************/
CavityPressure::~CavityPressure() {
}
