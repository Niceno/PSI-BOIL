#include "phasechange.h"

/******************************************************************************/
PhaseChange::PhaseChange(const Scalar & MDOT, 
                   const Scalar & TPR, 
                   const Scalar & TPRS,
                   const Scalar & CLR,
                   const Scalar & CLRS,
                   const Scalar & VS,
                   const Scalar & STEP,
                   const Vector & U, 
                   Times & T, 
                   Matter * f,
                   real r1, real r2,
                   Matter * s,
                   Nucleation * NUCL ) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( MDOT.domain(), MDOT, VS, & U, T, f, s, NULL ),
  tpr(&TPR),
  tprs(&TPRS),
  clr(&CLR),
  clrs(&CLRS),
  step(&STEP),
  dist  ( *MDOT.domain()),
  dflag ( *MDOT.domain()),
  iflag ( *MDOT.domain()),
  stmp  ( *MDOT.domain()),
  stmp2 ( *MDOT.domain()),
  delta( *MDOT.domain()),
  nx    ( *MDOT.domain()),
  ny    ( *MDOT.domain()),
  nz    ( *MDOT.domain()),
  txv   ( *MDOT.domain()),
  tyv   ( *MDOT.domain()),
  tzv   ( *MDOT.domain()),
  txl   ( *MDOT.domain()),
  tyl   ( *MDOT.domain()),
  tzl   ( *MDOT.domain()),
  M     ( *MDOT.domain())
{
  dist   = MDOT.shape();
  dflag  = MDOT.shape();
  iflag  = MDOT.shape();
  stmp   = MDOT.shape();
  stmp2  = MDOT.shape();
  delta = MDOT.shape();
  nx     = MDOT.shape();
  ny     = MDOT.shape();
  nz     = MDOT.shape();
  txv    = MDOT.shape();
  tyv    = MDOT.shape();
  tzv    = MDOT.shape();
  txl    = MDOT.shape();
  tyl    = MDOT.shape();
  tzl    = MDOT.shape();
  M      = MDOT.shape();
  nucl   = NUCL;

  /* set arguments */
  latent = r1;
  tsat = r2;
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
  phimin = 0.0;
  phimax = 1.0;
  phisurf = 0.5*(phimin+phimax);
  turbP = 0.9;

  epsl = 1.0e-1;
  //epsl = 1.0e-2;
  epsnorm = 1.0e-12;
  imodcal = 0;
  dxmin = dom->dxyz_min();

  /* set default values */
  Mmicro = 2.78E-05;
  Fmicro = -2.27E-02;
  use_int_res = false;
  resint = 0.0e0;

  alloc1d ( & hflux_total, 7);
  alloc1d ( & hflux_micro, 7);
  alloc1d ( & hflux_vapor, 7);
  alloc1d ( & area_sum, 7);
  alloc1d ( & area_l, 7);
  alloc1d ( & area_v, 7);

}	

/******************************************************************************/
PhaseChange::~PhaseChange() {
}
