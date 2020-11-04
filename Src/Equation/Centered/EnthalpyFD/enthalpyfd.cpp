#include "enthalpyfd.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
EnthalpyFD::EnthalpyFD(const Scalar & PHI, 
                       const Scalar & F,
                       const Vector & U,
                       const Vector & Uliq,
                       const Vector & Ugas,
                       Times & T,
                       Linear * S,
                       Matter * f,
                       const CommonHeatTransfer & CHT,
                       Matter * s) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S ),
  ftif   (  *PHI.domain()),
  cht(CHT),
  iflag(CHT.topo->iflag),
  iflagold(&(CHT.topo->iflagold)),
  uliq(&Uliq),
  ugas(&Ugas)
{
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());
  laminar=true;
  ao_conv = AccuracyOrder::Second();

  /* see header for explanation */
  if(solid()) {
    safe_solid = solid();
    accelerated_no_solid = false;
  } else {
    safe_solid = fluid();
    accelerated_no_solid = true;
  }

  ftif = phi.shape();
  phi.bnd_update();

  convection_set(TimeScheme::forward_euler());
  diffusion_set(TimeScheme::backward_euler());

  discretize();
}	
