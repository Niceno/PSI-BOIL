#include "nucleation.h"
#include "header.h"

/***************************************************************************//**
*  nucleation site model
*******************************************************************************/
Nucleation::Nucleation ( CommonHeatTransfer * CHT, Heaviside * HEAVI,
                         const Times * t, const real rs, 
                         Scalar * QSRC, const Sign SIG) :
  cht(CHT),
  heavi(HEAVI),
  matter_sig(SIG)
{
  vf = cht->topo->vf;
  clr= cht->topo->clr;
  qsrc=QSRC;
  time=t;
  rseed=rs;

  if(matter_sig==Sign::pos()) {
    rhol = cht->fluid()->rho(1);
    rhov = cht->fluid()->rho(0);
    lambdal = cht->fluid()->lambda(1);
    lambdav = cht->fluid()->lambda(0);
    cpl = cht->fluid()->cp(1);
    cpv = cht->fluid()->cp(0);
    mmass = cht->fluid()->mmass(0);
  } else {
    rhol = cht->fluid()->rho(0);
    rhov = cht->fluid()->rho(1);
    lambdal = cht->fluid()->lambda(0);
    lambdav = cht->fluid()->lambda(1);
    cpl = cht->fluid()->cp(0);
    cpv = cht->fluid()->cp(1);
    mmass = cht->fluid()->mmass(1);
  }
  latent = cht->fluid()->latent()->value();

#ifndef USE_VOF
  rcut = 4.*rseed;
#endif
  seed_period = 0.01;
  period_prevent_replant = 0.0;
  dxmin = vf->domain()->dxyz_min();
  eps = 1.5*dxmin;
  bzoning = false;
  zbtm = 0.0;
  threshold_c = 0.5;

  limit_zoning = true;
  zoning_limit_multiplier = 0.2;
}

/******************************************************************************/
Nucleation::~Nucleation() {
}
