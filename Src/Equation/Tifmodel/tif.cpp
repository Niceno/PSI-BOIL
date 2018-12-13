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
         const Scalar & MFLX,
         Matter * FLU,
         const Scalar * PRES,
         const Scalar * ADENS) :
  flu(FLU),
  mflx(&MFLX),  
  tif    (  *MFLX.domain()),
  tifold (  *MFLX.domain())
{

  pres = PRES;
  adens = ADENS;

  const int comp = 1; /* liquid */
  rhol = fluid()->rho(comp);
  tr = Tref;
  latent = Latent;
  mresis = Mresis;

  tif    = mflx.shape(); /* a mistake? */
  tifold = mflx.shape(); /* a mistake? */
  store_tif = false;
  variable_tif = true;
}
