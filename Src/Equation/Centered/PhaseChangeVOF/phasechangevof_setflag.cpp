#include "phasechangevof.h"
#include <cmath>

/******************************************************************************/
void PhaseChangeVOF::setflag() {
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
    if(clr[i][j][k]>=clrsurf){
      iflag[i][j][k]=3;
    } else {
      iflag[i][j][k]=-3;
    }
  }

#ifdef IB
  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      iflag[i][j][k]=-1001;
  }
#endif

#if 1
  for_ijk(i,j,k) {
    if(Interface(i,j,k)) {
      if(clr[i][j][k]<clrsurf)
        iflag[i][j][k]=-1;
      else
        iflag[i][j][k]=+1; /* clr = 0.5? questionable...*/
    }
  }
#else /* if you decide to use this one later on, check the enthalpy flagging
       * for a way how to efficiently treat the ibody.
       * Moreover, the interfacial flagging should really be unified among
       * classes...I am looking at you, Enthalpy, VOF, PhaseChange... */

  /* i-direction */
  for(int i=clr.si()-1; i<=clr.ei(); i++){
    for_vjk(clr,j,k){
       if((clr[i][j][k]-clrsurf)*(clr[i+1][j][k]-clrsurf)<=0.0){
         if(clr[i][j][k]<clrsurf){
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
      if((clr[i][j][k]-clrsurf)*(clr[i][j+1][k]-clrsurf)<=0.0){
        if(clr[i][j][k]<clrsurf){
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
       if((clr[i][j][k]-clrsurf)*(clr[i][j][k+1]-clrsurf)<=0.0){
        if(clr[i][j][k]<clrsurf){
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
#endif

  iflag.exchange_all();

  //boil::plot->plot(clr,iflag, "ib_setflag_clr-dflag", time->current_step());
  //exit(0);

  return;
}
