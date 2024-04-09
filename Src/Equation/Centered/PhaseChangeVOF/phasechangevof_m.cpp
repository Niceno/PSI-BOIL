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
      real qv = -lv*tnv[i][j][k];
      real ql =  ll*tnl[i][j][k];
      M[i][j][k] = (qv + ql) / latent;

      //boil::oout<<"PCV_m: "<<i<<" "<<j<<" "<<k<<" | "<<phi.zc(k)<<" | "<<clr[i][j][k]<<" "<<adens[i][j][k]<<" "<<fs[Comp::w()][i][j][k]<<" | "<<tnv[i][j][k]<<" "<<tnl[i][j][k]<<" | "<<M[i][j][k]<<boil::endl;
    } else {
      M[i][j][k] = 0.0;
    }
  }

  M.bnd_update();
  M.exchange_all();

  boil::timer.stop("phasechangevof m");

  return;
}

