#include "distance.h"

/******************************************************************************/
Distance::Distance(const Scalar & PHI, 
                   const Scalar & F,
                   const Vector & U, 
                   Times & T,
                   Krylov * S) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  v_fluid( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &v_fluid, NULL, S )  
{ 
  assert(PHI.domain() == F.domain());

  diffusion_set (TimeScheme::backward_euler());
  convection_set(TimeScheme::backward_euler());
}	

/******************************************************************************/
Distance::~Distance() {
}	

/*-----------------------------------------------------------------------------+
 '$Id: distance.cpp,v 1.2 2010/03/25 08:15:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
