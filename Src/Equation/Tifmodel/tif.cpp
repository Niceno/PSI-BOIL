#include "tif.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
TIF::TIF(const real Tref) {
  tr = Tref;
  variable_tif = false;
}

TIF::TIF(const real Tref, 
         const real Latent,
         const real Mresis,
         Matter * FLU,
         const Scalar & ADENS,
         const Scalar & MFLX,
         const Scalar * PRES) :
  flu(FLU),
  adens(&ADENS),  
  mflx(&MFLX),  
  tif    (*MFLX.domain()),
  tifold (*MFLX.domain())
{
  dpres = PRES;

  const int comp = 1; /* liquid */
  rhol = fluid()->rho(comp);
  tr = Tref;
  latent = Latent;
  mresis = Mresis;

  tmin = -boil::unreal;
  tmax = boil::unreal;
  weaklim = false;
  stronglim = false;
  clr = NULL;
  tpr = NULL;
  clrsurf = 0.5;

  tif    = mflx.shape(); /* a mistake? */
  tifold = mflx.shape(); /* a mistake? */

  store_tif = false;
  variable_tif = true;
  
  factor = 0.05;
  tint_field();
}
