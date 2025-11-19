#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::set_reconstruction_flag(const Scalar & scp, ScalarInt & flag,
                                        const int nlayer) {
/***************************************************************************//**
*  /brief crude flagging to set the reconstruction belt
*         output   :flag
*
*   flag:  -3 | -2 | -1 |  0 |  0 |  1 |  2 |  3 | ..   
*          ---+----+----+----+----+----+----+----+---
*   vapor                   |                        liquid
*                      free-surface
*
*   solid region: flag = -1001
*******************************************************************************/
  const int ifmax =  nlayer+4;
  const int ifmin = -nlayer-4;

  /*-----------+
  |  set flag  |
  +-----------*/
  for_aijk(i,j,k) {
    if(dom->ibody().off(i,j,k)) {
      flag[i][j][k]=-1001;
    } else {
      if(scp[i][j][k]<phisurf){
        flag[i][j][k]=ifmin;
      } else {
        flag[i][j][k]=ifmax;
      }
    }
  }

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=0; i<ni()-1; i++){
    for_jk(j,k){
      if ((scp[i][j][k]-phisurf)*(scp[i+1][j][k]-phisurf)<=0.0) {
        if(flag[i  ][j][k]>-1000) flag[i  ][j][k]=0;
        if(flag[i+1][j][k]>-1000) flag[i+1][j][k]=0;
      }
    }
  }
#if 0
  /* j-direction */
  for(int j=0; j<nj()-1; j++){
    for_ik(i,k){
      if((scp[i][j][k]-phisurf)*(scp[i][j+1][k]-phisurf)<=0.0){
        if(flag[i][j  ][k]>-1000) flag[i][j  ][k]=0;
        if(flag[i][j+1][k]>-1000) flag[i][j+1][k]=0;
      }
    }
  }
#endif
  /* k-direction */
  for(int k=0; k<nk()-1; k++){
    for_ij(i,j){
      if((scp[i][j][k]-phisurf)*(scp[i][j][k+1]-phisurf)<=0.0){
        if(flag[i][j][k  ]>-1000) flag[i][j][k  ]=0;
        if(flag[i][j][k+1]>-1000) flag[i][j][k+1]=0;
      }
    }
  }

  flag.exchange();

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
        if(abs(flag[i  ][j][k])==ifmax &&
           abs(flag[i+1][j][k])==(layer-1)){
           if(flag[i  ][j][k]>-1000)
             flag[i  ][j][k]=signum(layer,flag[i  ][j][k]);
        } else if(abs(flag[i+1][j][k])==ifmax &&
                  abs(flag[i  ][j][k])==(layer-1)){
           if(flag[i+1][j][k]>-1000)
             flag[i+1][j][k]=signum(layer,flag[i+1][j][k]);
        }
      }
    }
#if 0
    /* j-direction */
    for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
        if(abs(flag[i][j  ][k])==ifmax &&
           abs(flag[i][j+1][k])==(layer-1)){
           if(flag[i][j  ][k]>-1000)
             flag[i][j  ][k]=signum(layer,flag[i][j  ][k]);
        } else if(abs(flag[i][j+1][k])==ifmax &&
                  abs(flag[i][j  ][k])==(layer-1)){
           if(flag[i][j+1][k]>-1000)
             flag[i][j+1][k]=signum(layer,flag[i][j+1][k]);
        }
      }
    }
#endif
    /* k-direction */
    for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
        if(abs(flag[i][j][k  ])==ifmax &&
           abs(flag[i][j][k+1])==(layer-1)){
           if(flag[i][j][k  ]>-1000)
             flag[i][j][k  ]=signum(layer,flag[i][j][k  ]);
        } else if(abs(flag[i][j][k+1])==ifmax &&
                  abs(flag[i][j][k  ])==(layer-1)){
           if(flag[i][j][k+1]>-1000)
             flag[i][j][k+1]=signum(layer,flag[i][j][k+1]);
        }
      }
    }
    flag.exchange();
  }

#if 0
  boil::plot->plot(scp, flag, "scp-flag", time->current_step());
  //exit(0);
#endif

  return;
}
