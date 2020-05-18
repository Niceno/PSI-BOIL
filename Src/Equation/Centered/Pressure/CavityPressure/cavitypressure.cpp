#include "cavitypressure.h"

/******************************************************************************/
CavityPressure::CavityPressure(const Scalar & PHI,
                               const Scalar & F,
                               const Vector & U,
                               Times & T,
                               Linear * sm,
                               Matter * liq,
                               Topology & topo,
                               const real TENS,
                               const Scalar & CURV,
                               Sign SIG) :
  Pressure(PHI,F,U,T,sm,liq),
  fs(topo.fs),
  iflag(topo.iflag),
  sig(SIG),
  tens(TENS),
  clr(topo.clr),
  kappa(CURV)
{
  cavity_pressure = 0.0;

  clrsurf = 0.5;
  /* the in_gas function is initialized using a lambda expression */
  /* clr = 1 is liquid */
  if(SIG>0) {
    in_gas = [this](const real c) { return c<clrsurf; };
  /* clr = 0 is liquid */
  } else {
    in_gas = [this](const real c) { return c>clrsurf; }; 
  }

  /* to overwrite the default discretization imposed by parent constructor */
  discretize();

}	

/******************************************************************************/
CavityPressure::~CavityPressure() {
}
