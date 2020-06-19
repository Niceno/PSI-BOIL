#include "nucleation.h"

/***************************************************************************//**
*  nucleation site model
*******************************************************************************/
Nucleation::Nucleation ( Topology * TOPO, Heaviside * HEAVI,
                         const Scalar * TPR,
                         const Times * t, 
                         Matter * f, const real rs, 
                         Scalar * QSRC, const Sign SIG) :
  topo(TOPO),
  heavi(HEAVI),
  flu(f),
  matter_sig(SIG)
{
  vf = TOPO->vf;
  clr= TOPO->clr;
  tpr=TPR;
  qsrc=QSRC;
  time=t;
  rseed=rs;

  if(matter_sig==Sign::pos()) {
    rhol = fluid()->rho(1);
    rhov = fluid()->rho(0);
    lambdal = fluid()->lambda(1);
    lambdav = fluid()->lambda(0);
    cpl = fluid()->cp(1);
    cpv = fluid()->cp(0);
    mmass = fluid()->mmass(0);
  } else {
    rhol = fluid()->rho(0);
    rhov = fluid()->rho(1);
    lambdal = fluid()->lambda(0);
    lambdav = fluid()->lambda(1);
    cpl = fluid()->cp(0);
    cpv = fluid()->cp(1);
    mmass = fluid()->mmass(1);
  }
  latent = fluid()->latent()->value();

  rcut = 4.*rseed;
  seed_period = 0.01;
  period_cut_replant = 0.0001;
  dxmin = clr->domain()->dxyz_min();
  eps = 1.5*dxmin;
  bzoning = false;
  zbtm = 0.0;

  limit_zoning = true;
  zoning_limit_multiplier = 0.2;
}

/******************************************************************************/
Nucleation::~Nucleation() {
}
