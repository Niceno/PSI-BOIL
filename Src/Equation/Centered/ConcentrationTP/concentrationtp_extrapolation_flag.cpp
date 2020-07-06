#include "concentrationtp.h"

/******************************************************************************/
void ConcentrationTP::extrapolation_flag() {
/***************************************************************************//**
*  \brief set flag for eps extrapolation.
*         vapor cell  = -3
*         liquid cell = +3
*         int cell    = 0
*         vapor  cell next to int = -1
*         liquid cell next to int = +1
*         ibody  cell = -1000
*
*         ^ = bitwise XOR
*******************************************************************************/

  for_vijk(clr,i,j,k) {
    if((*topo->iflag)[i][j][k]>0){
      eflag[i][j][k]=3;
    } else {
      eflag[i][j][k]=-3;
    }
    if(heavi->status(i,j,k)==0) eflag[i][j][k] = 0;
  }

  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      eflag[i][j][k]=-1000;
  }

  eflag.bnd_update();
  eflag.exchange_all();

  //boil::plot->plot(clr,eflag, "clr-eflag", time->current_step());

  /* i-direction */
  for(int i=clr.si()-1; i<=clr.ei(); i++){
    for_vjk(clr,j,k){
      int & eflag_c  = eflag[i  ][j][k];
      int & eflag_p  = eflag[i+1][j][k];      
      bool bc = eflag_c;
      bool bp = eflag_p;
      if((bc^bp) && (eflag_c+eflag_p)>-500) {
        if(bc) {
          eflag_c = heavi->status(i  ,j,k);
        } else {
          eflag_p = heavi->status(i+1,j,k);
        }
      } 
    }
  }
  /* j-direction */
  for(int j=clr.sj()-1; j<=clr.ej(); j++){
    for_vik(clr,i,k){
      int & eflag_c  = eflag[i][j  ][k];
      int & eflag_p  = eflag[i][j+1][k]; 
      bool bc = eflag_c;
      bool bp = eflag_p;
      if((bc^bp) && (eflag_c+eflag_p)>-500) {
        if(bc) {
          eflag_c = heavi->status(i,j  ,k);
        } else {
          eflag_p = heavi->status(i,j+1,k);
        }  
      }
    }
  }
  /* k-direction */
  for(int k=clr.sk()-1; k<=clr.ek(); k++){
    for_vij(clr,i,j){
      int & eflag_c  = eflag[i][j][k  ];
      int & eflag_p  = eflag[i][j][k+1];
      bool bc = eflag_c;
      bool bp = eflag_p;
      if((bc^bp) && (eflag_c+eflag_p)>-500) {
        if(bc) {
          eflag_c = heavi->status(i,j,k  );
        } else {
          eflag_p = heavi->status(i,j,k+1);
        }
      }
    }
  }

  eflag.bnd_update();
  eflag.exchange_all();

#if 0
  if(time->current_step() > 0) {
    boil::plot->plot(clr,eflag, "clr-eflag", time->current_step());
    exit(0);
  }
#endif

  return;
}
