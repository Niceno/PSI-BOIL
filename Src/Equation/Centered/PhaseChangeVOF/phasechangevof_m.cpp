#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::m(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate M, usually in unit kg/m2s.
*         M = (qflux_liquid + qflux_vapor) / latent
*******************************************************************************/

  boil::timer.start("phasechangevof m");

#if 1
  for_ijk(i,j,k){
    if(Interface(i,j,k)){
      real lv = lambdav;
      real ll = lambdal;
      if (diff_eddy) {
        lv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
        ll += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real qv = -lv*tnv[i][j][k];
      real ql =  ll*tnl[i][j][k];
      M[i][j][k] = (qv + ql) / latent;
    } else {
      M[i][j][k] = 0.0;
    }
  }
#else
  
  real clr0 = 3.88968036e-13;
  real mflx0 = 0.063;

  real clrsum(0.0);
  real volsum(0.0);
  for_ijk(i,j,k) {
    clrsum += clr[i][j][k]*phi.dV(i,j,k);
    volsum += phi.dV(i,j,k);
  }
  boil::cart.sum_real(&clrsum);
  boil::cart.sum_real(&volsum);
  clrsum = (volsum - clrsum)/3.0/phi.dzc(1+boil::BW);
  clr0   = (volsum - clr0  )/3.0/phi.dzc(1+boil::BW);
  clrsum = sqrt(clrsum/boil::pi);
  clr0   = sqrt(clr0  /boil::pi);
  for_ijk(i,j,k) {
    if(Interface(i,j,k)){
      M[i][j][k] = mflx0*clr0/clrsum;
    } else {
      M[i][j][k] = 0.0;
    }
  }

#endif

#if 0
  real mflxval(0.0);
  real cnt(0.0);
  for_ijk(i,j,k) {
    if(Interface(i,j,k)){
      real locval = M[i][j][k];
      mflxval += locval;
      cnt += 1.0;
    }
  }
  boil::cart.sum_real(&mflxval);
  boil::cart.sum_real(&cnt);
  mflxval/=cnt;
  
  for_ijk(i,j,k){
    if(Interface(i,j,k)){
      M[i][j][k] = mflxval;
    } else {
      M[i][j][k] = 0.0;
    }
  }
#endif

  M.bnd_update();
  M.exchange_all();

  boil::timer.stop("phasechangevof m");

  return;
}

