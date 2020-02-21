#include "schrage.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
Schrage::Schrage(const real Tref, 
                 Matter * FLU,
                 const Scalar * ADENS,
                 const Scalar & MFLX,
                 const Scalar * PRES) :
  TIF(Tref,*ADENS),
  flu(FLU),
  mflx(&MFLX)
{
  adens = ADENS;
  dpres = PRES;

  rhol = fluid()->rho(1);
  real rhov = fluid()->rho(0);
  real mmass = fluid()->mmass(0);
  tr = Tref;
  mresis = pow(tr,1.5)/2.0/rhov/fluid()->latent()->value()
                      /sqrt(mmass/(2.0*boil::pi*boil::R));
  boil::oout<<"Schrage: Mass transfer resistance= "<<mresis<<boil::endl;

  variable_tif = true;
  
  tint_field();
}
