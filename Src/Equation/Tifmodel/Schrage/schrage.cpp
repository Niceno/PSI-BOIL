#include "schrage.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
Schrage::Schrage(const real Tref, 
                 Matter * FLU,
                 Topology * TOPO,
                 const Scalar & MFLX,
                 const Scalar * PRES) :
  TIF(Tref,TOPO),
  flu(FLU),
  mflx(&MFLX)
{
  dpres = PRES;

  rhol = fluid()->rho(1);
  real rhov = fluid()->rho(0);
  real mmass = fluid()->mmass(0);
  real latent = fluid()->latent()->value();
  tr = Tref;

  hresis = calculate_heat_transfer_resistance(tr,rhov,mmass,latent,1.);
  mresis = hresis * latent;
  boil::oout<<"Schrage: Transfer resistance= "<<hresis<<" "<<mresis<<boil::endl;

  variable_tif = true;

  boil::oout<<"Schrage: obsolete class. Exiting."<<boil::endl;
  exit(0);
  
  //tint_field();

}
/******************************************************************************/
/*----------------+
|  static members |
+----------------*/
real Schrage::calculate_heat_transfer_resistance(const real tr, const real rhov,
                                                 const real mmass,
                                                 const real latent,
                                                 const real accommodation) {
  const real accomult = accommodation;//2.*accommodation/(2.-accommodation);
  return std::pow(tr,1.5)/accomult/rhov/latent/latent
                         /sqrt(mmass/(2.0*boil::pi*boil::R));
}
