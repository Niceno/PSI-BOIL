#include "enthalpytif.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
EnthalpyTIF::EnthalpyTIF(const Scalar & PHI, 
                       const Scalar & F,
                       const Scalar & C,
                       const Vector & U,
                       Times & T,
                       Krylov * S,
                       Matter * f,
                       const real Tsat,
                       Matter * s,
                       const Vector * FS,
                       const real Latent,
                       const real Mresis,
                       const Scalar * MFLX,
                       const Scalar * PRES) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S ),
  clrold (  *C  .domain()),
  tif    (  *PHI.domain()),
  tifold (  *PHI.domain()),
  ftif   (  *PHI.domain()),
  ftifold(  *PHI.domain()),
  fdelta (  *PHI.domain()),
  iflag     (*C  .domain()),
  fsold(  *U  .domain())
{
  tsat = Tsat,
  rhol = fluid()->rho(1),
  rhov = fluid()->rho(0),
  cpl  = fluid()->cp(1),
  cpv  = fluid()->cp(0),
  lambdal = fluid()->lambda(1),
  lambdav = fluid()->lambda(0),
  clr = &C;
  clrsurf = 0.5;
  clrold = (*clr).shape();
  iflag  = (*clr).shape();
  store_clrold = false;
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());
  epsl=1.0e-2;
  turbP=0.9;
  laminar=true;

  fs = FS;

  if(fs) {
    for_m(m) {
      fsold(m) = (*fs)(m).shape();
    }
  }
 
  ftif = phi.shape();
  mflx = MFLX;
  pres = PRES;
  if(mflx||pres) {
    tif    = phi.shape();
    tifold = phi.shape();
    latent = Latent;
    mresis = Mresis;
    store_tif = false;
  }
  blendfactor = 0.05;

  phi.bnd_update();

  convection_set(TimeScheme::forward_euler());
  diffusion_set(TimeScheme::backward_euler());

  discretize();
}	

/******************************************************************************/
EnthalpyTIF::~EnthalpyTIF() {
}	

real EnthalpyTIF::blendfactor;

