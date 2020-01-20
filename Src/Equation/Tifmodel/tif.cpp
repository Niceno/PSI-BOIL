#include "tif.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
TIF::TIF(const real Tref) {
  tr = Tref;
  variable_tif = false;
}

TIF::TIF(const real Tref, 
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

  rhol = fluid()->rho(1);
  real rhov = fluid()->rho(0);
  real mmass = fluid()->mmass(0);
  tr = Tref;
  mresis = pow(tr,1.5)/2.0/rhov/fluid()->latent()->value()
                      /sqrt(mmass/(2.0*boil::pi*boil::R));
  boil::oout<<"TIFmodel: Mass transfer resistance= "<<mresis<<boil::endl;


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
