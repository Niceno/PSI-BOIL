#include "electpoten.h"

/******************************************************************************/
ElectPoten::ElectPoten(const Scalar & PHI, 
                   const Scalar & F, 
                   Vector & b,
                   Vector & j,
                   const Vector & U,
                   Times & T, 
                   Linear * sm,
                   Matter * f) :
  /* call parent's constructor. NULL is for solid */
  Centered( PHI.domain(), PHI, F, & U, T, f, NULL, sm ),
  B(&b),
  J(&j),
  uf( *PHI.domain()),
  vf( *PHI.domain()),
  wf( *PHI.domain())
{ 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == sm->domain());
  set_indices(B);
  set_indices(J);
  for_m(m) {
    uf(m) = (*u)(m).shape();
    vf(m) = (*u)(m).shape();
    wf(m) = (*u)(m).shape();
  }
  diffusion_set(TimeScheme::backward_euler());

  phi.bnd_update();
  discretize();
}	

/******************************************************************************/
ElectPoten::~ElectPoten() {
}
