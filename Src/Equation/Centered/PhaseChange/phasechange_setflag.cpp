#include "phasechange.h"
#include <cmath>

/******************************************************************************/
void PhaseChange::setflag() {
/***************************************************************************//**
*  \brief set flag.
          vapor cell  = -3
          nFSv cell   = -2
          FSv cell    = -1
          FSl cell    = 1
          nFSl cell   = 2
          liquid cell = 3
          ibody  cell = -1000
*******************************************************************************/

  for_vijk(clr,i,j,k) {
    if(clr[i][j][k]>=phisurf){
      iflag[i][j][k]=3;
    } else {
      iflag[i][j][k]=-3;
    }
  }
  /* i-direction */
  for(int i=clr.si()-1; i<=clr.ei(); i++){
    for_vjk(clr,j,k){
       if((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0){
         if(clr[i][j][k]<phisurf){
            iflag[i  ][j][k]=-1;
            iflag[i+1][j][k]=1;
         } else {
            iflag[i  ][j][k]=1;
            iflag[i+1][j][k]=-1;
         }
       }
    }
  }
  /* j-direction */
  for(int j=clr.sj()-1; j<=clr.ej(); j++){
    for_vik(clr,i,k){
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
        if(clr[i][j][k]<phisurf){
          iflag[i][j  ][k]=-1;
          iflag[i][j+1][k]=1;
        } else {
          iflag[i][j  ][k]=1;
          iflag[i][j+1][k]=-1;
        }
      }
    }
  }
  /* k-direction */
  for(int k=clr.sk()-1; k<=clr.ek(); k++){
    for_vij(clr,i,j){
       if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
        if(clr[i][j][k]<phisurf){
          iflag[i][j][k  ]=-1;
          iflag[i][j][k+1]=1;
        } else {
          iflag[i][j][k  ]=1;
          iflag[i][j][k+1]=-1;
        }
      }
    }
  }

  iflag.exchange_all();
  /* cells in vicinity of free surface cells */

  /* i-direction */
  for(int i=clr.si()-1; i<=clr.ei(); i++){
    for_vjk(clr,j,k){
      if(!(fabs(iflag[i][j][k])==1) != !(fabs(iflag[i+1][j][k])==1)){ /* xor */
        if(fabs(iflag[i][j][k]) == 1){
           iflag[i+1][j][k] = copysign(2,iflag[i  ][j][k]);
        } else {
           iflag[i  ][j][k] = copysign(2,iflag[i+1][j][k]);
        }
      }
    }
  }
  /* j-direction */
  for(int j=clr.sj()-1; j<=clr.ej(); j++){
    for_vik(clr,i,k){
      if(!(fabs(iflag[i][j][k])==1) != !(fabs(iflag[i][j+1][k])==1)){ /* xor */
        if(fabs(iflag[i][j][k]) == 1){
           iflag[i][j+1][k] = copysign(2,iflag[i][j  ][k]);
        } else {
           iflag[i][j  ][k] = copysign(2,iflag[i][j+1][k]);
        }
      }
    }
  }
  /* k-direction */
  for(int k=clr.sk()-1; k<=clr.ek(); k++){
    for_vij(clr,i,j){
      if(!(fabs(iflag[i][j][k])==1) != !(fabs(iflag[i][j][k+1])==1)){ /* xor */
        if(fabs(iflag[i][j][k]) == 1){
           iflag[i][j][k+1] = copysign(2,iflag[i][j][k  ]);
        } else {
           iflag[i][j][k  ] = copysign(2,iflag[i][j][k+1]);
        }
      }
    }
  }

#ifdef IB
  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      iflag[i][j][k]=-1001;
  }
#endif

  iflag.exchange_all();

  //boil::plot->plot(clr,iflag, "ib_setflag_clr-dflag", time->current_step());
  //exit(0);

  return;
}
