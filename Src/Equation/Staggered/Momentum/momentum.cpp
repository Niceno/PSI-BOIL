#include "momentum.h"

/******************************************************************************/
Momentum::Momentum(const Vector & U, 
                   const Vector & F,
                   Times & T,   
                   Linear * sm,
                   Matter * M) :
  Staggered(U.domain(), U, F, T, M, NULL, sm) { /* NULL is for solid */

  A[0] = new Matrix( u( Comp::u() ) );
  A[1] = new Matrix( u( Comp::v() ) );
  A[2] = new Matrix( u( Comp::w() ) );

  assert(U.domain() == F.domain());
  assert(U.domain() == sm->domain());

  /* find out which direction is pseudo */
  ifull = jfull = kfull = true;
  for_m(m)
    for(int b=0; b<u.bc(m).count(); b++) {
      if(u.bc(m).type(b) == BndType::pseudo()) {
        Dir d = u.bc(m).direction(b);
        if(d == Dir::imin() || d == Dir::imax())
          ifull = false;
        if(d == Dir::jmin() || d == Dir::jmax())
          jfull = false;
        if(d == Dir::kmin() || d == Dir::kmax())
          kfull = false;
      }
    }

  boil::oout<<"Momentum-full: "<<ifull<<" "<<jfull<<" "<<kfull<<boil::endl;
  ib_trust_vel_wall = false;

  u.bnd_update_nooutlet();

  discretize();

  v_phase_change=0.0;
}

/******************************************************************************/
Momentum::~Momentum() {
}
