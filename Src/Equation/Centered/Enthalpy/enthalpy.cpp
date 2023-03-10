#include "enthalpy.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
Enthalpy::Enthalpy(const Scalar & PHI, 
                   const Scalar & F,
                   const Vector & U, 
                   Times & T, 
                   Linear * S,
                   Matter * f,
                   Matter * s) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S )  
{ 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());

  turbPr = 0.9;  // default value of turbulent Prandtl number

  phi.bnd_update(); // must be called before discretization because of 
                    // dirichlet boundary condition etc.

  discretize();
}	

/******************************************************************************/
Enthalpy::~Enthalpy() {
}
