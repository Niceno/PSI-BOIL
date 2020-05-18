#include "phasechange4.h"

/******************************************************************************/
PhaseChange4::PhaseChange4(const Scalar & MDOT, 
                               const Scalar & MFLX,
                               const Scalar & TPR, 
                               const Scalar & TPRS,
                               const Scalar & VF,
                               const Scalar & VFS,
                               const Scalar & VS,
                               const Vector & U, 
                               Topology * TOPO,
                               const TIF & TIFMODEL,
                               Times & T, 
                               Matter * f,
                               Matter * s,
                               Sign SIG) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( MDOT.domain(), MDOT, VS, & U, T, f, s, NULL ),
  tpr(&TPR),
  tprs(&TPRS),
  clr(TOPO->clr),
  vf(&VF),
  vfs(&VFS),
  M(&MFLX),
  topo(TOPO),
  nx(TOPO->nx), /* these are aliases for easier use */
  ny(TOPO->ny),
  nz(TOPO->nz),
  adens(TOPO->adens),
  iflag(TOPO->iflag),
  bndtpr  ( *U   .domain() ),
  txv     ( *MDOT.domain()),
  tyv     ( *MDOT.domain()),
  tzv     ( *MDOT.domain()),
  txl     ( *MDOT.domain()),
  tyl     ( *MDOT.domain()),
  tzl     ( *MDOT.domain()),
  tnl     ( *MDOT.domain()),
  tnv     ( *MDOT.domain()),
  tempflag( *MDOT.domain()),
  tifmodel(TIFMODEL),
  sig(SIG)
{
  txv      = MDOT.shape();
  tyv      = MDOT.shape();
  tzv      = MDOT.shape();
  txl      = MDOT.shape();
  tyl      = MDOT.shape();
  tzl      = MDOT.shape();
  tnl      = MDOT.shape();
  tnv      = MDOT.shape();
  tempflag = MDOT.shape();
  for_m(m) bndtpr(m) = U(m).shape(); /* a mistake? */

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

  /* set arguments */
  rhol = fluid()->rho(1);
  rhov = fluid()->rho(0);
  lambdal = fluid()->lambda(1);
  lambdav = fluid()->lambda(0);
  cpl = fluid()->cp(1);
  cpv = fluid()->cp(0);

  /* set constants */
  clrsurf = 0.5;
  turbP = 0.9;

  /* this should be consistent with enthalpy class */
  //epsl = 5.0e-2;
  epsl = 1.0e-2;

  /* flags */
  use_second_order_accuracy = true;
  use_unconditional_extrapolation = false;
  discard_points_near_interface = true;

}	

/******************************************************************************/
PhaseChange4::~PhaseChange4() {
}	

