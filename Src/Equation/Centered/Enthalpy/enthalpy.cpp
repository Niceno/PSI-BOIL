#include "enthalpy.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
Enthalpy::Enthalpy(const Scalar & PHI, 
                   const Scalar & F,
                   const Vector & U, 
                   Times & T, 
                   Krylov * S,
                   Matter * f,
                   Matter * s) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S )  
{ 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());

  phi.bnd_update(); // must be called before discretization because of 
                    // dirichlet boundary condition etc.

  discretize();
}	

/******************************************************************************/
Enthalpy::~Enthalpy() {
}
