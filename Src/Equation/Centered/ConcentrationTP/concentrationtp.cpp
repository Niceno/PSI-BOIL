#include "concentrationtp.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
ConcentrationTP::ConcentrationTP(const Scalar & PHI,
                                 const Scalar & F,
                                 const Vector & U,
                                 const Vector & FLUXCLR, 
                                 Heaviside * HEAVI,
                                 Topology * TOPO,
                                 Times & T, 
                                 Linear * S,
                                 Matter * f,     /* the diffusing species */
                                 const Sign SIG
                                ) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, NULL, S ), /* NULL is for solid */
  /* color is actually volume fraction */
  clr(TOPO->vf),
  clrold(&(TOPO->vfold)),
  eflag ( *PHI.domain()),
  colorflow(&FLUXCLR),
  heavi(HEAVI),
  topo(TOPO),
  sig(SIG)
{
  rho_dif = (f->rho());     /* pointer at property */
  dcoef   = (f->gamma());   /* pointer at property */

  eflag   = phi.shape();
 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());
  turbS = 0.7;
  laminar = true;

  phi.bnd_update();
  /* must be called before discretization because of
     dirichlet boundary condition etc. */
  convection_set(TimeScheme::forward_euler());
  diffusion_set(TimeScheme::backward_euler());

  discretize();
}	

/******************************************************************************/
ConcentrationTP::~ConcentrationTP() {
}	
