#include "phasechangevof.h"

/******************************************************************************/
PhaseChangeVOF::PhaseChangeVOF(const Scalar & MDOT, 
                               const Scalar & TPR, 
                               const Scalar & TPRS,
                               const Scalar & CLR,
                               const Scalar & CLRS,
                               const Scalar & VS,
                               const Vector & U, 
                               const Scalar & NX,
                               const Scalar & NY,
                               const Scalar & NZ,
                               const Vector & FS,
                               Times & T, 
                               Matter * f,
                               real LAT, real TS,
                               Scalar * TIF,
                               Matter * s) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( MDOT.domain(), MDOT, VS, & U, T, f, s, NULL ),
  tpr(&TPR),
  tprs(&TPRS),
  clr(&CLR),
  clrs(&CLRS),
  mx(&NX),
  my(&NY),
  mz(&NZ),
  fs(&FS),
  nx(*NX.domain()),
  ny(*NY.domain()),
  nz(*NZ.domain()),
  bndtpr ( *U   .domain() ),
  M      ( *MDOT.domain()),
  txv    ( *MDOT.domain()),
  tyv    ( *MDOT.domain()),
  tzv    ( *MDOT.domain()),
  txl    ( *MDOT.domain()),
  tyl    ( *MDOT.domain()),
  tzl    ( *MDOT.domain()),
  stmp   ( *MDOT.domain()),
  stmp2  ( *MDOT.domain()),
  delta  ( *MDOT.domain()),
  iflag  ( *MDOT.domain()),
  gradclr( *CLR .domain())
{
  M       = MDOT.shape();
  txv     = MDOT.shape();
  tyv     = MDOT.shape();
  tzv     = MDOT.shape();
  txl     = MDOT.shape();
  tyl     = MDOT.shape();
  tzl     = MDOT.shape();
  stmp    = MDOT.shape();
  stmp2   = MDOT.shape();
  delta   = MDOT.shape();
  iflag   = MDOT.shape();
  gradclr = MDOT.shape();
  nx  = NX.shape();
  ny  = NY.shape();
  nz  = NZ.shape();
  for_m(m) bndtpr(m) = U(m).shape(); /* a mistake? */


  tif = TIF;

  /* set arguments */
  latent = LAT;
  tsat = TS;
  rhol = fluid()->rho(1);
  rhov = fluid()->rho(0);
  rhoave = 0.5*(rhol+rhov);
  rhodlt = fabs(rhol-rhov);
  lambdal = fluid()->lambda(1);
  lambdav = fluid()->lambda(0);
  cpl = fluid()->cp(1);
  cpv = fluid()->cp(0);

  /* set constants */
  tempnull = -1000.0;
  pi = acos(-1.0);
  phisurf = 0.5;
  turbP = 0.9;

  epsl = 1.0e-1;
  //epsl = 1.0e-2;
  dxmin = dom->dxyz_min();

#if 1
  upwind_flag = false;
#else
  upwind_flag = true;
#endif
}	

/******************************************************************************/
PhaseChangeVOF::~PhaseChangeVOF() {
}	

