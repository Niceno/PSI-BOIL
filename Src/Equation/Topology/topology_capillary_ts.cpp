#include "topology.h"

/* static version, also a static member */
real Topology::capillary_ts(const Matter & mixed,
                            const real dx, const real coef)  {
  return coef*sqrt(dx*dx*dx*(mixed.rho(1)+mixed.rho(0))
                  /mixed.sigma()->value());
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
  for_vijk(vfold,i,j,k) {
    if(interface(i,j,k)) {
      real uval = 0.5*(vel[Comp::u()][i][j][k]+vel[Comp::u()][i+1][j][k]);
      real vval = 0.5*(vel[Comp::v()][i][j][k]+vel[Comp::v()][i][j+1][k]);
      real wval = 0.5*(vel[Comp::w()][i][j][k]+vel[Comp::w()][i][j][k+1]);
      real un = uval*(*nx)[i][j][k]+vval*(*ny)[i][j][k]+wval*(*nz)[i][j][k];
      real usq = uval*uval+vval*vval+wval*wval;
      real ut = std::max(usq-un*un,0.0); /* this should be guaranteed... */
      ut = sqrt(ut);
      ut = ut/std::min(vfold.dxc(i),std::min(vfold.dyc(j),vfold.dzc(k)));
      if(ut>utmax)
        utmax = ut;
    }
  }
  boil::cart.max_real(&utmax);

  return 1./(1./capillary_ts(mixed,dxmin,coef)+utmax);
}

/* phasic velocities */
/******************************************************************************/
real Topology::capillary_ts(const Matter & mixed, const Vector & vel1, 
                            const Vector & vel2,  const real coef) const {
/******************************************************************************/

  real utmax(0.);
  for_vijk(vfold,i,j,k) {
    if(interface(i,j,k)) {
      real uval1 = 0.5*(vel1[Comp::u()][i][j][k]+vel1[Comp::u()][i+1][j][k]);
      real vval1 = 0.5*(vel1[Comp::v()][i][j][k]+vel1[Comp::v()][i][j+1][k]);
      real wval1 = 0.5*(vel1[Comp::w()][i][j][k]+vel1[Comp::w()][i][j][k+1]);
      real uval2 = 0.5*(vel2[Comp::u()][i][j][k]+vel2[Comp::u()][i+1][j][k]);
      real vval2 = 0.5*(vel2[Comp::v()][i][j][k]+vel2[Comp::v()][i][j+1][k]);
      real wval2 = 0.5*(vel2[Comp::w()][i][j][k]+vel2[Comp::w()][i][j][k+1]);
      real un1 = uval1*(*nx)[i][j][k]+vval1*(*ny)[i][j][k]+wval1*(*nz)[i][j][k];
      real un2 = uval2*(*nx)[i][j][k]+vval2*(*ny)[i][j][k]+wval2*(*nz)[i][j][k];
      real usq1 = uval1*uval1+vval1*vval1+wval1*wval1;
      real usq2 = uval2*uval2+vval2*vval2+wval2*wval2;
      real ut1 = std::max(usq1-un1*un1,0.0); /* this should be guaranteed... */
      real ut2 = std::max(usq2-un2*un2,0.0); 
      real ut = sqrt(std::max(ut1,ut2));
      ut = ut/std::min(vfold.dxc(i),std::min(vfold.dyc(j),vfold.dzc(k)));
      if(ut>utmax)
        utmax = ut;
    }
  }
  boil::cart.max_real(&utmax);

  return 1./(1./capillary_ts(mixed,dxmin,coef)+utmax);
}

