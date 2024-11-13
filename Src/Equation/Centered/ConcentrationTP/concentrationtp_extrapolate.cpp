#include "concentrationtp.h"

/******************************************************************************/
void ConcentrationTP::extrapolate() {
/***************************************************************************//**
*  \brief Extrapolate eps across the interface.
*******************************************************************************/

  /* flagging (different from topology flag) */
  extrapolation_flag();

  /* extrapolate epsilon from vapor to liquid */
#if 0
  topo->extrapolate(phi,-matter_sig,{-matter_sig},eflag);
#else
  for_avijk(phi,i,j,k) {
    if(matter_sig*eflag[i][j][k]<0) {
      // full-of-liquid cells
      phi[i][j][k]=0.0;
      eflag2[i][j][k] = 0;  // eps will be extrapolated if eflag2 = 0
    } else {
      // interface-cells and full-of-vapor cells
      eflag2[i][j][k] = 1;  // eps will be fixed if eflag2 = 1
    }
  }

  eflag = eflag2;
  stmp = phi;

  #if 0
    /* visualize flag */
    for_ijk(i,j,k){
      stmp[i][j][k]=eflag[i][j][k];
    }
    boil::plot->plot(clr,phi,stmp, "before-c-eps-eflag", time->current_step());
    //exit(0);
    stmp = phi;
  #endif

  for(int iloop=1; iloop<3; iloop++) { 
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) continue;
      if(eflag[i][j][k]==0) {
        int inb = std::min(1,eflag[i-1][j][k]) + std::min(1,eflag[i+1][j][k])
                + std::min(1,eflag[i][j-1][k]) + std::min(1,eflag[i][j+1][k])
                + std::min(1,eflag[i][j][k-1]) + std::min(1,eflag[i][j][k+1]);
        if(inb >= 1) {
          stmp[i][j][k] = (real(std::min(1,eflag[i-1][j][k])) * phi[i-1][j][k]
                        +  real(std::min(1,eflag[i+1][j][k])) * phi[i+1][j][k]
                        +  real(std::min(1,eflag[i][j-1][k])) * phi[i][j-1][k]
                        +  real(std::min(1,eflag[i][j+1][k])) * phi[i][j+1][k]
                        +  real(std::min(1,eflag[i][j][k-1])) * phi[i][j][k-1]
                        +  real(std::min(1,eflag[i][j][k+1])) * phi[i][j][k+1])
                        /real(inb);
          if(stmp[i][j][k]<0.0||stmp[i][j][k]>1.0) {
            std::cout<<"ConcentrationTP::extrapolate::eps= "<<stmp[i][j][k]<<"\n";
            std::cout<<"at proc,i,j,k="<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
            std::cout<<"inb= "<<std::min(1,eflag[i-1][j][k])<<" "<<std::min(1,eflag[i+1][j][k])<<" "
                <<std::min(1,eflag[i][j-1][k]) + std::min(1,eflag[i][j+1][k])<<" "
                <<std::min(1,eflag[i][j][k-1])<<" "<<  std::min(1,eflag[i][j][k+1])<<"\n";
            std::cout<<"eps= "<<phi[i-1][j][k]<<" "<<phi[i+1][j][k]<<" "<<phi[i][j-1][k]<<" "
                              <<phi[i][j+1][k]<<" "<<phi[i][j][k-1]<<" "<<phi[i][j][k+1]<<"\n";
            std::cout<<"clr= "<<clr[i-1][j][k]<<" "<<clr[i+1][j][k]<<" "<<clr[i][j-1][k]<<" "
                              <<clr[i][j+1][k]<<" "<<clr[i][j][k-1]<<" "<<clr[i][j][k+1]<<"\n";
            exit(0);
          }
          eflag2[i][j][k] = 2;  /* eflag=2 for extrapolated */
        }
      }
    }
    stmp.bnd_update(); 
    eflag2.bnd_update();
    stmp.exchange();
    eflag2.exchange();
    phi = stmp;
    eflag = eflag2;
  }

  #if 0
    /* visualize flag */
    for_ijk(i,j,k){
      stmp[i][j][k]=eflag[i][j][k];
    }
    boil::plot->plot(clr,phi,stmp, "after-c-eps-eflag", time->current_step());
    exit(0);
  #endif
#endif

  // limit eps between 0 and 1
  for_ijk(i,j,k){
    phi[i][j][k]=std::min(1.0,std::max(0.0,phi[i][j][k]));
  }

  return;
}
