#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::m(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate M, usually in unit kg/m2s.
*         M = (qflux_liquid + qflux_vapor) / latent
*******************************************************************************/

  boil::timer.start("phasechangevof m");

  for_ijk(i,j,k){
    if(Interface(i,j,k)){
      real lv = lambdav;
      real ll = lambdal;
      if (diff_eddy) {
        lv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
        ll += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
#if 0
      real qv = -lv * ( txv[i][j][k]*nx[i][j][k]
                      + tyv[i][j][k]*ny[i][j][k]
                      + tzv[i][j][k]*nz[i][j][k]);
      real ql =  ll * ( txl[i][j][k]*nx[i][j][k]
                      + tyl[i][j][k]*ny[i][j][k]
                      + tzl[i][j][k]*nz[i][j][k]);
#else 
      real qv = -lv*tnv[i][j][k];
      real ql =  ll*tnl[i][j][k];
#endif
      M[i][j][k] = (qv + ql) / latent;
    } else {
      M[i][j][k] = 0.0;
    }
  }

#if 1
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

  M.exchange_all();

  boil::timer.stop("phasechangevof m");

  return;
}

