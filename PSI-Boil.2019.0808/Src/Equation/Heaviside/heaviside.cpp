#include "heaviside.h"

Heaviside::Heaviside(const Scalar * CLR, Scalar * PHI, Scalar * ADENS,
                     const real CLRSURF) : 
  clr(CLR), dom((*CLR).domain()), phi(PHI), adens(ADENS),
  flag(*CLR->domain()), clrsurf(CLRSURF), bflag_struct(*CLR)
{
  flag.copy_shape(clr->shape());

  for( int b=0; b<clr->bc().count(); b++ ) {
    if(    clr->bc().type(b) == BndType::dirichlet()
        || clr->bc().type(b) == BndType::inlet()
        || clr->bc().type(b) == BndType::outlet()
        || clr->bc().type(b) == BndType::insert()
        || clr->bc().type(b) == BndType::convective()
       ) {
      flag.bc().type(b) = BndType::neumann();
      boil::oout << "Adjusting b.c.s for flag at " << b
                << boil::endl;
    }
  }

  /* for vertex interpolation reasons */
  assert(boil::nano>boil::pico);
}
