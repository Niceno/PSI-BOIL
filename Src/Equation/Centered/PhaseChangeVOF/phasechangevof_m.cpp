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
      real qv = -lv * ( txv[i][j][k]*nx[i][j][k]
                      + tyv[i][j][k]*ny[i][j][k]
                      + tzv[i][j][k]*nz[i][j][k]);
      real ql =  ll * ( txl[i][j][k]*nx[i][j][k]
                      + tyl[i][j][k]*ny[i][j][k]
                      + tzl[i][j][k]*nz[i][j][k]);
      M[i][j][k] = (qv + ql) / latent;
    } else {
      M[i][j][k] = 0.0;
    }
  }
  M.exchange_all();

  boil::timer.stop("phasechangevof m");

  return;
}

