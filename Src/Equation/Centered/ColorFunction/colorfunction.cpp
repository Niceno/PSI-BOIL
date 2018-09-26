#include "colorfunction.h"

/******************************************************************************/
ColorFunction::ColorFunction(const Scalar & PHI, 
                             const Scalar & F,
                             const real & con, 
                             const real & den, 
                             const Vector & U, 
                             Times & T,
                             Krylov * S) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  v_fluid( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &v_fluid, NULL, S ), 
  nx( *PHI.domain() ),
  ny( *PHI.domain() ),
  nz( *PHI.domain() ),
  ndiv( *PHI.domain() ),
  surf( *PHI.domain() ),
  phi_old( *PHI.domain() )
/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  nx      = phi.shape();
  ny      = phi.shape();
  nz      = phi.shape();
  ndiv    = phi.shape();
  surf    = phi.shape();
  phi_old = phi.shape();

  assert(PHI.domain() == F.domain());

  diffusion_set (TimeScheme::backward_euler());
  //convection_set(TimeScheme::backward_euler());
  convection_set(TimeScheme::adams_bashforth());

  discretize();
}	

/******************************************************************************/
ColorFunction::~ColorFunction() {
}	

/*-----------------------------------------------------------------------------+
 '$Id: colorfunction.cpp,v 1.27 2014/08/06 08:37:03 sato Exp $'/
+-----------------------------------------------------------------------------*/
