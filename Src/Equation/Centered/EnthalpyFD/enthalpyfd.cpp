#include "enthalpyfd.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
EnthalpyFD::EnthalpyFD(const Scalar & PHI, 
                       const Scalar & F,
                       const Scalar & C,
                       const Vector & U,
                       Times & T,
                       Krylov * S,
                       Matter * f,
                       TIF & tintmodel,
                       Matter * s,
                       const Vector * FS,
                       const Scalar * ADENS) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S ),
  clrold (  *C  .domain()),
  ftif   (  *PHI.domain()),
  ftifold(  *PHI.domain()),
  adensold( *PHI.domain()),
  iflag     (*C  .domain()),
  fsold(  *U  .domain()),
  tifmodel(tintmodel)
{
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
  adens = ADENS;

  if(fs) {
    for_m(m) {
      fsold(m) = (*fs)(m).shape();
    }
  }
 
  ftif = phi.shape();
  if(adens)
    adensold = (*adens).shape();

  phi.bnd_update();

  convection_set(TimeScheme::forward_euler());
  diffusion_set(TimeScheme::backward_euler());

  discretize();
}	

/******************************************************************************/
EnthalpyFD::~EnthalpyFD() {
}	

/******************************************************************************/
