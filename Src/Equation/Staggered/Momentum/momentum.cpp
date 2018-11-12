#include "momentum.h"

/******************************************************************************/
Momentum::Momentum(const Vector & U, 
                   const Vector & F,
                   Times & T,   
                   Krylov * sm,
                   Matter * M) :
  Staggered(U.domain(), U, F, T, M, NULL, sm) { /* NULL is for solid */

  A[0] = new Matrix( u( Comp::u() ) );
  A[1] = new Matrix( u( Comp::v() ) );
  A[2] = new Matrix( u( Comp::w() ) );

  assert(U.domain() == F.domain());
  assert(U.domain() == sm->domain());

  insert_bc();

  discretize();

  v_phase_change=0.0;
}

/******************************************************************************/
Momentum::~Momentum() {
}
