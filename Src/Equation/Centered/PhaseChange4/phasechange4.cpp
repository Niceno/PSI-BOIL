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
                               HTWallModel * HTM,
                               Sign SIG) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( MDOT.domain(), MDOT, VS, & U, T, f, s, NULL ),
  tpr(&TPR),
  tprs(&TPRS),
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
  htwallmodel(HTM),
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

  /* heat transfer wall model should be equal to the one of enthalpy,
   * that's why it is a pointer. however, the default argument is a
   * nullptr, we need to fix that */
  if(!htwallmodel) {
    default_value_for_htwallmodel = true;
    htwallmodel = new HTWallModel();
  } else {
    default_value_for_htwallmodel = false;
  }

  /* set arguments */
  rhol = fluid()->rho(1);
  rhov = fluid()->rho(0);
  lambdal = fluid()->lambda(1);
  lambdav = fluid()->lambda(0);
  cpl = fluid()->cp(1);
  cpv = fluid()->cp(0);

  /* set constants */
  turbP = 0.9;

  /* flags */
  accuracy_order = 2;
  use_unconditional_extrapolation = false;

  /* false-true: stability/accuracy tradeoff */
  discard_points_near_interface = false; //true;

}	

/******************************************************************************/
PhaseChange4::~PhaseChange4() {
  if(default_value_for_htwallmodel)
    delete htwallmodel;
}	

