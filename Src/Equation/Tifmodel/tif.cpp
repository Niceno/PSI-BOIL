#include "tif.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
TIF::TIF(const real Tref) {
  tr = Tref;
  variable_tif = false;

  /* unused within constant model */
  tmin = -boil::unreal;
  tmax = boil::unreal;
  weaklim = false;
  store_tif = false;
  factor = 0.05;
}

/* initialises tif scalars */
TIF::TIF(const real Tref, Topology * TOPO) :
  tif(*TOPO->adens->domain()),
  tifold(*TOPO->adens->domain()),
  topo(TOPO),
  iflag(TOPO->iflag),
  tempflag(*TOPO->adens->domain()),
  tempflag2(*TOPO->adens->domain()),
  stmp(*TOPO->adens->domain())
{
  tr = Tref;
  variable_tif = false;

  tmin = -boil::unreal;
  tmax = boil::unreal;
  weaklim = false;
  store_tif = false;
  factor = 0.05;

  tif    = (*TOPO->adens).shape(); 
  tifold = (*TOPO->adens).shape();
  tempflag  = (*TOPO->adens).shape();
  tempflag2 = (*TOPO->adens).shape();
  stmp = (*TOPO->adens).shape();
}

/******************************************************************************/
/*----------------+
|  static members |
+----------------*/
real TIF::calculate_heat_transfer_resistance(const real tr, const real rhov,
                                             const real mmass,
                                             const real latent,
                                             const real accommodation) {
  const real accomult = accommodation;//2.*accommodation/(2.-accommodation);
  return std::pow(tr,1.5)/accomult/rhov/latent/latent
                         /sqrt(mmass/(2.0*boil::pi*boil::R));
}

