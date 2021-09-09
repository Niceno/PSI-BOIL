#include "antoine.h"

/***************************************************************************//**
*  Constructors
*******************************************************************************/
Antoine::Antoine(const real Tref, 
                 Topology * TOPO,
                 const Scalar & EPS,
                 const real a,
                 const real b,
                 const real c) :
  TIF(Tref,TOPO),
  eps(&EPS)
{
  A = a;
  B = b;
  C = c;
  /* conversion to Kelvin */
  C_K = C-273.15;
  tr = Tref;

  variable_tif = true;
  
  //tint_field();
}
