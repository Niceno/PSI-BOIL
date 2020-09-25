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
                       Topology * TOPO,
                       TIF & tintmodel,
                       Matter * s,
                       HTWallModel * HTM) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, s, S ),
  ftif   (  *PHI.domain()),
  topo(TOPO),
  iflag(TOPO->iflag),
  iflagold(&(TOPO->iflagold)),
  uliq(&Uliq),
  ugas(&Ugas),
  tifmodel(tintmodel),
  htwallmodel(HTM)
{
  rhol = fluid()->rho(1),
  rhov = fluid()->rho(0),
  cpl  = fluid()->cp(1),
  cpv  = fluid()->cp(0),
  lambdal = fluid()->lambda(1),
  lambdav = fluid()->lambda(0),
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());
  epsl=1.0e-2; /* this appears in diff_matrix but should not play a role */
  turbP=0.9;
  laminar=true;

  /* see header for explanation */
  if(solid()) {
    safe_solid = solid();
    accelerated_no_solid = false;
  } else {
    safe_solid = fluid();
    accelerated_no_solid = true;
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

  ftif = phi.shape();
  phi.bnd_update();

  convection_set(TimeScheme::forward_euler());
  diffusion_set(TimeScheme::backward_euler());

  discretize();
}	

/******************************************************************************/
EnthalpyFD::~EnthalpyFD() {
  if(default_value_for_htwallmodel)
    delete htwallmodel;
}	

/******************************************************************************/
