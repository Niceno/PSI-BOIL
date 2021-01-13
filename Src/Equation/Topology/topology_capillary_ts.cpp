#include "topology.h"

real Topology::wave_vel(const Matter & mixed,
                        const real dx, const real coef) const {
  return sqrt(mixed.sigma()->value()
            /((mixed.rho(1)+mixed.rho(0))*dx))/coef;
}

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
  
  //real inv_stat = 1./capillary_ts(mixed,dxmin,coef);
  real umax(0.);
  for_vijk(vfold,i,j,k) {
    if(interface(i,j,k)) {
      real uval = 0.5*(vel[Comp::u()][i][j][k]+vel[Comp::u()][i+1][j][k]);
      real vval = 0.5*(vel[Comp::v()][i][j][k]+vel[Comp::v()][i][j+1][k]);
      real wval = 0.5*(vel[Comp::w()][i][j][k]+vel[Comp::w()][i][j][k+1]);
      real un = uval*(*nx)[i][j][k]+vval*(*ny)[i][j][k]+wval*(*nz)[i][j][k];
      real usq = uval*uval+vval*vval+wval*wval;
      real ut = std::max(usq-un*un,0.0); /* this should be guaranteed... */
      ut = sqrt(ut);
      real dxloc = std::min(vfold.dxc(i),std::min(vfold.dyc(j),vfold.dzc(k)));
#if 1
      ut = (ut+wave_vel(mixed,dxloc))/dxloc;
      if(ut>umax)
        umax = ut;
#else
  #if 0
      ut = ut/dxloc;
      un = un/dxloc;
      ut += inv_stat;
      usq = ut*ut + un*un;
      usq = sqrt(usq);
  #else
      ut += wave_vel(mixed,dxloc);
      usq = ut*ut + un*un;
      usq = sqrt(usq)/dxloc;
  #endif
      if(usq>umax)
        umax = usq;
#endif 
    }
  }
  boil::cart.max_real(&umax);

  //return 1./(inv_stat+umax);
  return 1./umax;
}

/* phasic velocities */
/******************************************************************************/
real Topology::capillary_ts(const Matter & mixed, const Vector & vel1, 
                            const Vector & vel2,  const real coef) const {
/******************************************************************************/

  //real inv_stat = 1./capillary_ts(mixed,dxmin,coef);
  real umax(0.);
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
      ut1 = sqrt(ut1);
      ut2 = sqrt(ut2);
      real dxloc = std::min(vfold.dxc(i),std::min(vfold.dyc(j),vfold.dzc(k)));
#if 1
      real ut = std::max(ut1,ut2); 
      ut = (ut+wave_vel(mixed,dxloc))/dxloc;
      if(ut>umax)
        umax = ut;
#else
  #if 0
      ut1 = ut1/dxloc;
      un1 = un1/dxloc;
      ut2 = ut2/dxloc;
      un2 = un2/dxloc;
      ut1 += inv_stat;
      ut2 += inv_stat;
      usq1 = ut1*ut1 + un1*un1;
      usq2 = ut2*ut2 + un2*un2;
  #else
      ut1 += wave_vel(mixed,dxloc);
      ut2 += wave_vel(mixed,dxloc);
      usq1 = ut1*ut1 + un1*un1;
      usq2 = ut2*ut2 + un2*un2;
      usq1 = sqrt(usq1)/dxloc;
      usq2 = sqrt(usq2)/dxloc;
  #endif
      real usq = std::max(usq1,usq2);
      if(usq>umax)
        umax = usq;
#endif
    }
  }
  boil::cart.max_real(&umax);

  //return 1./(inv_stat+umax);
  return 1./umax;
}

