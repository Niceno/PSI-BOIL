#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::m() {
/***************************************************************************//**
*  \brief calculate M, usually in unit kg/m2s.
*         M = (qflux_liquid + qflux_vapor) / latent
*******************************************************************************/

  for_ijk(i,j,k) {
    if(interface(i,j,k)) {
      real qv = -tnv[i][j][k];
      real ql =  tnl[i][j][k];
      M[i][j][k] = (qv + ql) / fluid()->latent(i,j,k);

#if 0
      if(i==3&&j==3&&k==3) {
      boil::aout<<"PC4_m: "<<i<<" "<<j<<" "<<k<<" | "<<phi.xc(i)<<" | "<<clr[i][j][k]<<" "<<adens[i][j][k]<<" "<<(*(topo->fs))[Comp::k()][i][j][k]<<" "<<(*(topo->fs))[Comp::k()][i+1][j][k]<<" | "<<tpr[i][j][k]<<" "<<tpr[i][j][k-1]<<" "<<tpr[i][j][k+1]<<" | "<<tpr[i-1][j][k]<<" "<<tpr[i+1][j][k]<<" | "<<qv<<" "<<ql<<" "<<qv+ql<<" | "<<tzv[i][j][k]<<" "<<nx[i][j][k]<<" "<<nz[i][j][k]<<boil::endl;
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

