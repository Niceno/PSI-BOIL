#include "nucleation.h"

/***************************************************************************//**
*  nucleation site model
*******************************************************************************/
Nucleation::Nucleation ( Scalar * c, Scalar *tp, Scalar *qs,
                     const Times * t, Scalar & DM,
                     Matter *f, const real rs, const real dm, const real l,
                     const real ca, const real tst, Matter *s ) :
  dmicro(&DM),
  flu(f),
  sol(s)
{
  clr=c;
  tpr=tp;
  qsrc=qs;
  time=t;
  rseed=rs;
  dmicro_min = dm;
  latent = l;
  cang = ca;
  tsat = tst;

  /* default value */
  dmicro = boil::exa;
  seed_period = 0.001;
  dxmin       = c->domain()->dxyz_min();
  zbtm        = 0.0;
  store_dSprev = false;
  slope = 4.46e-3;  // Utaka's coefficient for water
  boil::oout<<"Nucleation:slope= "<<slope<<"\n";
  b_slope = 0.0;
  exp_slope = 1.0;
  rmax = 1.0e+300;
  bzoning = false;
  boptdat = false;
  eps_clr = 1.0e-4;
  clrsurf = 0.5;
  range_zoning = 1.0;
  set_heat_sink(true);
  set_pre_heat_sink(false);
  threshold_c = 0.5;

  /* material properties */
  rhov =    fluid()->rho(0);
  lambdav = fluid()->lambda(0);
  cpv =     fluid()->cp(0);
  cps =     solid()->cp()->value();

  /* allocate */
  const int n = boil::maxi(clr->ni(),clr->nj(),clr->nk());
  dSprev = new real*[n];
  for (int i=0; i<n; i++) {
    dSprev[i] = new real[n];
  }

}

/******************************************************************************/
Nucleation::~Nucleation() {
}
