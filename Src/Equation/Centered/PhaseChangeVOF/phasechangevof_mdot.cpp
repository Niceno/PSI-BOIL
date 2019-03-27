#include "phasechangevof.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChangeVOF::mdot() {
/***************************************************************************//**
*         output : phi [kg/m3s]
*******************************************************************************/

#if 0
  boil::oout<<"PCVOF::mdot"<<boil::endl;
#endif

  for_ijk(i,j,k){
   if(Interface(i,j,k)){
      if(dom->ibody().on(i,j,k)){
        real mdotc=M[i][j][k];

        /* iso-surface area */
        phi[i][j][k] = mdotc * adens[i][j][k];
        //phi[i][j][k] = mdot_cut(phi[i][j][k],clr[i][j][k]);
      } else {
        phi[i][j][k] = 0.0;
      }
    } else {
      phi[i][j][k] = 0.0;
    }
  }
#if 1
  phi.exchange_all();
#else
  /* cut */
  mdot_cut();
#endif

#if 0
  for_vijk(phi,i,j,k) {
    if(phi[i][j][k]>17300.&&phi[i][j][k]<17400.) {
 #if 0
    
      real qv = -( txv[i][j][k]*nx[i][j][k]
                      + tyv[i][j][k]*ny[i][j][k]
                      + tzv[i][j][k]*nz[i][j][k]);
      real ql =  ( txl[i][j][k]*nx[i][j][k]
                      + tyl[i][j][k]*ny[i][j][k]
                      + tzl[i][j][k]*nz[i][j][k]);
      boil::aout<<"PCV::mdot "<<i<<" "<<j<<" "<<k<<" | "<<phi[i][j][k]
                           <<" | "<<txv[i][j][k]<<" "<<txl[i][j][k]
                           <<" | "<<tyv[i][j][k]<<" "<<tyl[i][j][k]
                           <<" | "<<tzv[i][j][k]<<" "<<tzl[i][j][k]
                           <<" | "<<iflag[i][j][k]<<" "<<qv<<" "<<ql
                           <<" | "<<nx[i][j][k]<<" "<<ny[i][j][k]<<" "<<nz[i][j][k]
                           <<boil::endl;
    }
  #else
      real tm0;
      boil::aout<<"PCV::mdot "<<i<<" "<<j<<" "<<k<<" | "<<phi[i][j][k]
                           <<" | "<<txl[i][j][k]
                           <<" | "<<tpr[i-1][j][k]<<" "<<tpr[i][j][k]<<" "<<tpr[i+1][j][k]
                           //<<" | "<<gradtz9(+1,i,j,k)<<" "<<gradtz8(+1,i,j,k)
                           //<<" | "<<bndtpr[Comp::w()][i][j][k]<<" "<<bndtpr[Comp::w()][i][j][k+1]
                           //<<" | "<<gradtz9(-1,i,j,k)<<" "<<gradtz9(+1,i,j,k)
                           //<<" | "<<gradtz8(-1,i,j,k)<<" "<<gradtz8(+1,i,j,k)
                           //<<" | "<<distance_z(i,j,k,-1,tm0)<<" "<<distance_z(i,j,k,+1,tm0)
                           //<<" | "<<distance_z(i,j,k,-1,tm0)/clr.dzc(k)<<" "<<distance_z(i,j,k,+1,tm0)/clr.dzc(k)
                           <<boil::endl;

    }
  #endif
  }
#endif

  return;
}

