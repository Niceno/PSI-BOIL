#include "phasefield.h"

/******************************************************************************/
PhaseField::PhaseField(const Scalar & PHI, 
                       const Scalar & F,
                       const real & con, 
                       const real & den, 
                       const Vector & U, 
                       Times & T,
                       Krylov * S) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  jelly( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &jelly, NULL, S ), 
  phi_old( *PHI.domain() )
/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  phi_old       = phi.shape();

  assert(PHI.domain() == F.domain());

  /* needed?
  diffusion_set (TimeScheme::backward_euler());
  convection_set(TimeScheme::backward_euler());
  */
}	

/******************************************************************************/
PhaseField::~PhaseField() {
}	

/*-----------------------------------------------------------------------------+
 '$Id: phasefield.cpp,v 1.6 2011/03/29 13:40:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
