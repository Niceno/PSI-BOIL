#include "concentrationtp.h"

/******************************************************************************/
void ConcentrationTP::extrapolate() {
/***************************************************************************//**
*  \brief Extrapolate eps across the interface.
*******************************************************************************/

  /* flagging (different from topology flag) */
  extrapolation_flag();

  /* extrapolate epsilon from vapor to liquid */
#if 1
  //boil::oout<<"extrapolate: "<<phi[14][33][10]<<"\n";
  topo->extrapolate(phi,-matter_sig,{-matter_sig},eflag);
  //topo->extrapolate(phi,matter_sig,{-matter_sig},eflag);
  //topo->extrapolate(phi,-matter_sig,{matter_sig},eflag);
  //topo->extrapolate(phi,matter_sig,{matter_sig},eflag);
  //boil::oout<<"extrapolate: "<<phi[14][33][10]<<"\n";
#else
  for_avijk(phi,i,j,k) {
    if(matter_sig*eflag[i][j][k]<0) {
      phi[i][j][k]=0.0;
      eflag2[i][j][k] = 0;
    } else {
      eflag2[i][j][k] = 1;
    }
  }
  eflag = eflag2;

  stmp = phi;

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
  if(time->current_step()==1) {
    /* visualize flag */
    for_ijk(i,j,k){
      stmp[i][j][k]=eflag[i][j][k];
    }
    boil::plot->plot(clr,phi,stmp, "c-eps-eflag", time->current_step());
    exit(0);
  }
  #endif
#endif

  return;
}
