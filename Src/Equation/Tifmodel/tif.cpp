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
TIF::TIF(const real Tref, const Topology & topo) :
  tif(*topo.adens->domain()),
  tifold(*topo.adens->domain()),
  iflag(topo.iflag),
  tempflag(*topo.adens->domain()),
  tempflag2(*topo.adens->domain()),
  stmp(*topo.adens->domain())
{
  tr = Tref;
  variable_tif = false;

  tmin = -boil::unreal;
  tmax = boil::unreal;
  weaklim = false;
  store_tif = false;
  factor = 0.05;

  tif    = (*topo.adens).shape(); 
  tifold = (*topo.adens).shape();
  tempflag  = (*topo.adens).shape();
  tempflag2 = (*topo.adens).shape();
  stmp = (*topo.adens).shape();
}
