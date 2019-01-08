#include "phasechangevof.h"

/******************************************************************************/
PhaseChangeVOF::PhaseChangeVOF(const Scalar & MDOT, 
                               const Scalar & MFLX,
                               const Scalar & TPR, 
                               const Scalar & TPRS,
                               const Scalar & CLR,
                               const Scalar & CLRS,
                               const Scalar & VS,
                               const Vector & U, 
#if 0
                               const Scalar & NX,
                               const Scalar & NY,
                               const Scalar & NZ,
                               const Scalar & ADENS,
                               const Vector & FS,
#else
                               const VOF & vof,
#endif
                               const TIF & TIFMODEL,
                               Times & T, 
                               Matter * f,
                               real LAT,
                               Matter * s) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( MDOT.domain(), MDOT, VS, & U, T, f, s, NULL ),
  tpr(&TPR),
  tprs(&TPRS),
  clr(&CLR),
  clrs(&CLRS),
  M(&MFLX),
#if 0
  mx(&NX),
  my(&NY),
  mz(&NZ),
  fs(&FS),
  adens(&ADENS),
  nx(*NX.domain()),
  ny(*NY.domain()),
  nz(*NZ.domain()),
#else
  mx(&(vof.nx)),
  my(&(vof.ny)),
  mz(&(vof.nz)),
  fs(&(vof.fs)),
  adens(&(vof.adens)),
  nx(*(vof.nx).domain()),
  ny(*(vof.ny).domain()),
  nz(*(vof.nz).domain()),
#endif
  bndtpr ( *U   .domain() ),
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
  tifmodel(TIFMODEL)
{
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
#if 0
  nx  = NX.shape();
  ny  = NY.shape();
  nz  = NZ.shape();
#else
  nx  = (vof.nx).shape();
  ny  = (vof.ny).shape();
  nz  = (vof.nz).shape();
#endif
  for_m(m) bndtpr(m) = U(m).shape(); /* a mistake? */

  /* set arguments */
  latent = LAT;
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

  epsl = 5.0e-2;
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

