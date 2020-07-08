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
                                 const real critvf,
                                 const bool heavivolume,
                                 const Sign SIG
                                ) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, NULL, S ), /* NULL is for solid */
  /* color is actually volume fraction */
  clr(TOPO->vf),
  clrold(&(TOPO->vfold)),
  vfold ( *PHI.domain()),
  stmp  ( *PHI.domain()),
  eflag ( *PHI.domain()),
  eflag2( *PHI.domain()),
  colorflow(&FLUXCLR),
  heavi(HEAVI),
  topo(TOPO),
  col_crit(critvf),
  use_heaviside_instead_of_vf(heavivolume),
  matter_sig(SIG)
{
  rho_dif = (f->rho());     /* pointer at property */
  dcoef   = (f->gamma());   /* pointer at property */

  vfold   = phi.shape();
  stmp    = phi.shape();
  //eflag   = phi.shape(); /* don't use, BndCnd would be the same pointer! */
  //eflag2  = phi.shape();
  eflag.copy_shape(phi.shape());
  eflag2.copy_shape(phi.shape());
 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());
  turbS = 0.7;
  laminar = true;
  store_vf = false;

  for( int b=0; b<phi.bc().count(); b++ ) {
    if(    phi.bc().type(b) == BndType::dirichlet()
        || phi.bc().type(b) == BndType::inlet()
        || phi.bc().type(b) == BndType::outlet()
        || phi.bc().type(b) == BndType::insert()
        || phi.bc().type(b) == BndType::convective()
       ) {
       eflag.bc().type(b) = BndType::neumann();
       eflag2.bc().type(b) = BndType::neumann();
       boil::oout << "Adjusting b.c.s for flags at " << b
                  << boil::endl;
    }
  }

  /* volume fraction or heaviside */
  if(use_heaviside_instead_of_vf) {
    vfval = [this](const int i, const int j, const int k)
                  { return heavi->vf(i,j,k); };
    vfvalold = [this](const int i, const int j, const int k)
                     { return vfold[i][j][k]; };
  } else {
    vfval = [this](const int i, const int j, const int k)
                  { return std::min(1.0,std::max(0.0,clr[i][j][k])); };
    vfvalold = [this](const int i, const int j, const int k)
                     { return std::min(1.0,std::max(0.0,clrold[i][j][k])); };
  }

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
