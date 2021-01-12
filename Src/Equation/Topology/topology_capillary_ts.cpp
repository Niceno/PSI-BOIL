#include "topology.h"

/* static version */
real Topology::capillary_ts(const Matter & mixed, const real coef) const {
  return coef*sqrt(dxmin*dxmin*dxmin*(mixed.rho(1)+mixed.rho(0))
                  /mixed.sigma());
}

/******************************************************************************/
real Topology::capillary_ts(const Matter & mixed, const Vector & vel, 
                            const real coef) const {
/***************************************************************************//**
*  \brief Calculate dynamic time capillary time step according to
*         Denner & van Wachem,JCP,285(2015),24-40
*         Results: ts_limit
*******************************************************************************/

  real utmax(0.);
  for_ijk(i,j,k) {
    if(Interface(i,j,k)) {
      real uval = 0.5*(vel[Comp::u()][i][j][k]+vel[Comp::u()][i+1][j][k]);
      real vval = 0.5*(vel[Comp::v()][i][j][k]+vel[Comp::v()][i][j+1][k]);
      real wval = 0.5*(vel[Comp::w()][i][j][k]+vel[Comp::w()][i][j][k+1]);
      real un = uval*mx[i][j][k]+vval*my[i][j][k]+wval*mz[i][j][k];
      real usq = uval*uval+vval*vval+wval*wval;
      real ut = std::max(usq-un*un,0.0); /* this should be guaranteed... */
      ut = sqrt(ut);
      ut = ut/std::min(phi.dxc(i),std::min(phi.dyc(j),phi.dzc(k)));
      if(ut>utmax)
        utmax = ut;
    }
  }
  boil::cart.max_real(&utmax);

  return 1./(1./capillary_ts(mixed,coef)+utmax);
}

/* phasic velocities */
/******************************************************************************/
real Topology::capillary_ts(const Matter & mixed, const Vector & vel1, 
                            const Vector & vel2,  const real coef) const {
/******************************************************************************/

  real utmax(0.);
  for_ijk(i,j,k) {
    if(Interface(i,j,k)) {
      real uval1 = 0.5*(vel1[Comp::u()][i][j][k]+vel1[Comp::u()][i+1][j][k]);
      real vval1 = 0.5*(vel1[Comp::v()][i][j][k]+vel1[Comp::v()][i][j+1][k]);
      real wval1 = 0.5*(vel1[Comp::w()][i][j][k]+vel1[Comp::w()][i][j][k+1]);
      real uval2 = 0.5*(vel2[Comp::u()][i][j][k]+vel2[Comp::u()][i+1][j][k]);
      real vval2 = 0.5*(vel2[Comp::v()][i][j][k]+vel2[Comp::v()][i][j+1][k]);
      real wval2 = 0.5*(vel2[Comp::w()][i][j][k]+vel2[Comp::w()][i][j][k+1]);
      real un1 = uval1*mx[i][j][k]+vval1*my[i][j][k]+wval1*mz[i][j][k];
      real un2 = uval2*mx[i][j][k]+vval2*my[i][j][k]+wval2*mz[i][j][k];
      real usq1 = uval1*uval1+vval1*vval1+wval1*wval1;
      real usq2 = uval2*uval2+vval2*vval2+wval2*wval2;
      real ut1 = std::max(usq1-un1*un1,0.0); /* this should be guaranteed... */
      real ut2 = std::max(usq2-un2*un2,0.0); 
      real ut = sqrt(std::max(ut1,ut2));
      ut = ut/std::min(phi.dxc(i),std::min(phi.dyc(j),phi.dzc(k)));
      if(ut>utmax)
        utmax = ut;
    }
  }
  boil::cart.max_real(&utmax);

  return 1./(1./capillary_ts(mixed,coef)+utmax);
}

