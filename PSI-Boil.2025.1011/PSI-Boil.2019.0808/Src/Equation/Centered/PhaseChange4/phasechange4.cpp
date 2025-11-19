#include "phasechange4.h"

/******************************************************************************/
PhaseChange4::PhaseChange4(const Scalar & MDOT, 
                           const Scalar & MFLX,
                           const Scalar & TPRS,
                           const Scalar & VFS,
                           const Scalar & VS,
                           const Vector & U, 
                           const CommonHeatTransfer & CHT,
                           Times & T, 
                           Matter * f,
                           Matter * s,
                           Sign SIG) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( MDOT.domain(), MDOT, VS, & U, T, f, s, NULL ),
  tprs(&TPRS),
  vfs(&VFS),
  M(&MFLX),
  cht(CHT),
  vf(CHT.topo->vf), /* these are aliases for easier use */
  nx(CHT.topo->nx),
  ny(CHT.topo->ny),
  nz(CHT.topo->nz),
  txv     ( *MDOT.domain()),
  tyv     ( *MDOT.domain()),
  tzv     ( *MDOT.domain()),
  txl     ( *MDOT.domain()),
  tyl     ( *MDOT.domain()),
  tzl     ( *MDOT.domain()),
  tnl     ( *MDOT.domain()),
  tnv     ( *MDOT.domain()),
  matter_sig(SIG)
{
#if 0 /* don't use this, it creates BndCnd pointers */
  txv      = MDOT.shape();
  tyv      = MDOT.shape();
  tzv      = MDOT.shape();
  txl      = MDOT.shape();
  tyl      = MDOT.shape();
  tzl      = MDOT.shape();
  tnl      = MDOT.shape();
  tnv      = MDOT.shape();
#else
  txv.copy_shape(MDOT.shape());
  tyv.copy_shape(MDOT.shape());
  tzv.copy_shape(MDOT.shape());
  txl.copy_shape(MDOT.shape());
  tyl.copy_shape(MDOT.shape());
  tzl.copy_shape(MDOT.shape());
  tnl.copy_shape(MDOT.shape());
  tnv.copy_shape(MDOT.shape());
#endif

  for( int b=0; b<phi.bc().count(); b++ ) {
    if(    phi.bc().type(b) == BndType::dirichlet()
        || phi.bc().type(b) == BndType::inlet()
        || phi.bc().type(b) == BndType::insert()
        || phi.bc().type(b) == BndType::convective()
       ) {
       txv.bc().type(b) = BndType::neumann();
       tyv.bc().type(b) = BndType::neumann();
       tzv.bc().type(b) = BndType::neumann();
       txl.bc().type(b) = BndType::neumann();
       tyl.bc().type(b) = BndType::neumann();
       tzl.bc().type(b) = BndType::neumann();

       boil::oout << "Adjusting b.c.s for temperature gradients at " << b
                  << boil::endl;
    }
  }

  /* flags */
  accuracy_order = AccuracyOrder::Second();
  use_unconditional_extrapolation = false;

  /* false-true: stability/accuracy tradeoff */
  discard_points_near_interface = false;

}	

/******************************************************************************/
PhaseChange4::~PhaseChange4() {
}	

