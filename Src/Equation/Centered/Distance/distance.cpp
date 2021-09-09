#include "distance.h"

/******************************************************************************/
Distance::Distance(const Scalar & PHI, 
                   const Scalar & F,
                   const Vector & U, 
                   Times & T,
                   Linear * S) :
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
