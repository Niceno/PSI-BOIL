#include "pressure.h"

/******************************************************************************/
Pressure::Pressure(const Scalar & PHI, 
                   const Scalar & F, 
                   const Vector & U,
                   Times & T, 
                   Krylov * sm,
                   Matter * f) :
  /* call parent's constructor. NULL is for solid */
  Centered( PHI.domain(), PHI, F, & U, T, f, NULL, sm ) 
{ 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == sm->domain());

  diffusion_set(TimeScheme::backward_euler());

  discretize();
}	

/******************************************************************************/
Pressure::~Pressure() {
}

/*-----------------------------------------------------------------------------+
 '$Id: pressure.cpp,v 1.32 2011/03/29 13:40:26 niceno Exp $'/
+-----------------------------------------------------------------------------*/
