#include "phasechangevof.h"

/******************************************************************************/
PhaseChangeVOF::PhaseChangeVOF(const Scalar & MDOT, 
                               const Scalar & MFLX,
                               const Scalar & TPR, 
                               const Scalar & TPRS,
                               const Scalar & VF,
                               const Scalar & VFS,
                               const Scalar & VS,
                               const Vector & U, 
                               Topology & topo,
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
  clr(topo.clr),
  vf(&VF),
  vfs(&VFS),
  M(&MFLX),
  nx(topo.nx),
  ny(topo.ny),
  nz(topo.nz),
  fs(topo.fs),
  adens(topo.adens),
  iflag(topo.iflag),
  bndtpr  ( *U   .domain() ),
  txv     ( *MDOT.domain()),
  tyv     ( *MDOT.domain()),
  tzv     ( *MDOT.domain()),
  txl     ( *MDOT.domain()),
  tyl     ( *MDOT.domain()),
  tzl     ( *MDOT.domain()),
  tnl     ( *MDOT.domain()),
  tnv     ( *MDOT.domain()),
  stmp    ( *MDOT.domain()),
  stmp2   ( *MDOT.domain()),
  delta   ( *MDOT.domain()),
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
  stmp     = MDOT.shape();
  stmp2    = MDOT.shape();
  delta    = MDOT.shape();
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
  rhoave = 0.5*(rhol+rhov);
  rhodlt = fabs(rhol-rhov);
  lambdal = fluid()->lambda(1);
  lambdav = fluid()->lambda(0);
  cpl = fluid()->cp(1);
  cpv = fluid()->cp(0);

  /* set constants */
  pi = acos(-1.0);
  clrsurf = 0.5;
  turbP = 0.9;

  tol_ext = 1e-7;
  epsl = 5.0e-2;
  //epsl = 1.0e-2;
  dxmin = dom->dxyz_min();

  near_wall_modelling = false;

#if 1
  upwind_flag = false;
#else
  upwind_flag = true;
#endif
}	

/******************************************************************************/
PhaseChangeVOF::~PhaseChangeVOF() {
}	

