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
TIF::TIF(const real Tref, const Scalar & s) :
  tif(*s.domain()),
  tifold(*s.domain()) {
  tr = Tref;
  variable_tif = false;

  /* unused within constant model */
  tmin = -boil::unreal;
  tmax = boil::unreal;
  weaklim = false;
  store_tif = false;
  factor = 0.05;

  tif    = s.shape(); 
  tifold = s.shape();
}
