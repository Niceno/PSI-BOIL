#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::m(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate M, usually in unit kg/m2s.
*         M = (qflux_liquid + qflux_vapor) / latent
*******************************************************************************/

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
      M[i][j][k] = (qv + ql) / fluid()->latent(i,j,k);

#if 0
      if(i==3&&j==3&&k==3) {
      boil::aout<<"PCV_m: "<<i<<" "<<j<<" "<<k<<" | "<<phi.xc(i)<<" | "<<clr[i][j][k]<<" "<<adens[i][j][k]<<" "<<fs[Comp::k()][i][j][k]<<" "<<fs[Comp::k()][i+1][j][k]<<" | "<<tpr[i][j][k]<<" "<<tpr[i][j][k-1]<<" "<<tpr[i][j][k+1]<<" | "<<tpr[i-1][j][k]<<" "<<tpr[i+1][j][k]<<" | "<<qv<<" "<<ql<<" "<<qv+ql<<" | "<<tzv[i][j][k]<<" "<<nx[i][j][k]<<" "<<nz[i][j][k]<<boil::endl;
      }
#endif
    } else {
      M[i][j][k] = 0.0;
    }
  }

  M.bnd_update();
  M.exchange_all();

  return;
}

