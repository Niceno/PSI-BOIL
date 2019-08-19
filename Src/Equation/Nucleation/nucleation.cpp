#include "nucleation.h"

/***************************************************************************//**
*  nucleation site model
*******************************************************************************/
Nucleation::Nucleation ( Scalar * c, Scalar *tp, Scalar *qs,
                     const Times * t, Scalar & DM,
                     Matter *f, const real rs, const real dm, const real l,
                     const real ca ) :
  dmicro(&DM),
  flu(f)
{
  clr=c;
  tpr=tp;
  qsrc=qs;
  time=t;
  rseed=rs;
  dmicro_min = dm;
  latent = l;
  cang = ca;

  /* default value */
  dmicro = boil::exa;
  seed_period = 0.01;
  period_cut_replant = 0.0001;
  dxmin       = c->domain()->dxyz_min();
  zbtm        = 0.0;
  store_dSprev = false;
  slope = 4.46e-3;  // Utaka's coefficient for water
  exp_slope = 1.0;
  boil::oout<<"Nucleation:slope= "<<slope<<"\n";
  rmax = 1.0e+300;
  bzoning = false;

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
