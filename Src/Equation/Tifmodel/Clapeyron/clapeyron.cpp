#include "clapeyron.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
Clapeyron::Clapeyron(const real Tref, 
                     const Scalar * ADENS,
                     const Scalar & EPS,
                     /* the following should be changed later on to matter */
                     const real MM,
                     const real LATENT,
                     const real LATENT_SLP) :
  TIF(Tref,*ADENS),
  eps(&EPS)
{
  adens = ADENS;

  latent_cst = LATENT;
  latent_slp = LATENT_SLP;
  Rm = boil::R/MM;
  tr = Tref;
  tri = 1./Tref;

  latent = latent_cst + latent_slp * tr;

  errmax = 1e-3;
  nmax = 4;

  variable_tif = true;
  
  tint_field();
}