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
*******************************************************************************/

  for_avijk((*clr),i,j,k) {
    if((*topo->iflag)[i][j][k]>0){
      eflag[i][j][k]=3;
    } else {
      eflag[i][j][k]=-3;
    }
    if(heavi->status(i,j,k)==0) eflag[i][j][k] = 0;
  }

  for_aijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      eflag[i][j][k]=-1000;
  }

  //boil::plot->plot((*clr),eflag, "clr-eflag", time->current_step());

  /* i-direction */
  for(int i=(*clr).si()-1; i<=(*clr).ei(); i++){
    for_vjk((*clr),j,k){
      bool eflagc  = eflag[i  ][j][k];
      bool eflagp  = eflag[i+1][j][k];      
      if(eflagc ^ eflagp && (eflagc+eflagp)>-500) {
        if(eflagc) {
          eflag[i  ][j][k] = heavi->status(i  ,j,k); /* should be zero or one */
        } else {
          eflag[i+1][j][k] = heavi->status(i+1,j,k); /* should be zero or one */
        }
      } 
    }
  }
  /* j-direction */
  for(int j=(*clr).sj()-1; j<=(*clr).ej(); j++){
    for_vik((*clr),i,k){
      bool eflagc  = eflag[i][j  ][k];
      bool eflagp  = eflag[i][j+1][k]; 
      if(eflagc ^ eflagp && (eflagc+eflagp)>-500) {
        if(eflagc) {
          eflag[i][j  ][k] = heavi->status(i,j  ,k); /* should be zero or one */
        } else {
          eflag[i][j+1][k] = heavi->status(i,j+1,k); /* should be zero or one */
        }  
      }
    }
  }
  /* k-direction */
  for(int k=(*clr).sk()-1; k<=(*clr).ek(); k++){
    for_vij((*clr),i,j){
      bool eflagc  = eflag[i][j][k  ];
      bool eflagp  = eflag[i][j][k+1];
      if(eflagc ^ eflagp && (eflagc+eflagp)>-500) {
        if(eflagc) {
          eflag[i][j][k  ] = heavi->status(i,j,k  ); /* should be zero or one */
        } else {
          eflag[i][j][k+1] = heavi->status(i,j,k+1); /* should be zero or one */
        }
      }
    }
  }

  eflag.exchange_all();

#if 0
  if(time->current_step() > 0) {
    boil::plot->plot((*clr),eflag, "clr-eflag", time->current_step());
    exit(0);
  }
#endif

  return;
}
