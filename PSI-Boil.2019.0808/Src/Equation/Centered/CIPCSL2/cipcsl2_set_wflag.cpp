#include "cipcsl2.h"
#include <iomanip>
using namespace std;

/******************************************************************************/
void CIPCSL2::set_wflag() {
/***************************************************************************//**
*  /brief flagging
*           output   :wflag
*
*   wflag: -3 | -2 | -1 |  0 |  0 |  1 |  2 |  3 | ..   
*          ---+----+----+----+----+----+----+----+---
*   vapor                   |                        liquid
*                      free-surface
*
*   solid region       : wflag = -1001
*   wall adjacent cells: wflag = -1001
*******************************************************************************/
  const int ifmax =  nlayer+4;
  const int ifmin = -nlayer-4;

  /*-----------+
  |  set flag  |
  +-----------*/
  for_aijk(i,j,k) {
    if(clr[i][j][k]<phisurf){
      wflag[i][j][k]=ifmin;
    } else {
      wflag[i][j][k]=ifmax;
    }
  }

  /* solid */
  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      wflag[i][j][k]=-1001;
  }

#ifdef IB
  /* wall adjacent cells */
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
      wflag[i][j][k]=-1001;
  }
#endif

  insert_bc_flag(wflag, true);
  wflag.exchange();

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=0; i<ni()-1; i++){
    for_jk(j,k){
      if ((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0) {
        if(wflag[i  ][j][k]>-1000) wflag[i  ][j][k]=0;
        if(wflag[i+1][j][k]>-1000) wflag[i+1][j][k]=0;
      }
    }
  }
  /* j-direction */
  for(int j=0; j<nj()-1; j++){
    for_ik(i,k){
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
        if(wflag[i][j  ][k]>-1000) wflag[i][j  ][k]=0;
        if(wflag[i][j+1][k]>-1000) wflag[i][j+1][k]=0;
      }
    }
  }
  /* k-direction */
  for(int k=0; k<nk()-1; k++){
    for_ij(i,j){
      if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
        if(wflag[i][j][k  ]>-1000) wflag[i][j][k  ]=0;
        if(wflag[i][j][k+1]>-1000) wflag[i][j][k+1]=0;
      }
    }
  }
  insert_bc_flag(wflag, true);
  wflag.exchange();

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
        if(int(abs(wflag[i  ][j][k]))==ifmax &&
           int(abs(wflag[i+1][j][k]))==(layer-1)){
           if(wflag[i  ][j][k]>-1000)
             wflag[i  ][j][k]=layer*int(copysign(1.0,wflag[i  ][j][k]));
        } else if(int(abs(wflag[i+1][j][k]))==ifmax &&
                  int(abs(wflag[i  ][j][k]))==(layer-1)){
           if(wflag[i+1][j][k]>-1000)
             wflag[i+1][j][k]=layer*int(copysign(1.0,wflag[i+1][j][k]));
        }
      }
    }
    /* j-direction */
    for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
        if(int(abs(wflag[i][j  ][k]))==ifmax &&
           int(abs(wflag[i][j+1][k]))==(layer-1)){
           if(wflag[i][j  ][k]>-1000)
             wflag[i][j  ][k]=layer*int(copysign(1.0,wflag[i  ][j][k]));
        } else if(int(abs(wflag[i][j+1][k]))==ifmax &&
                  int(abs(wflag[i][j  ][k]))==(layer-1)){
           if(wflag[i][j+1][k]>-1000)
             wflag[i][j+1][k]=layer*int(copysign(1.0,wflag[i][j+1][k]));
        }
      }
    }
    /* k-direction */
    for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
        if(int(abs(wflag[i][j][k  ]))==ifmax &&
           int(abs(wflag[i][j][k+1]))==(layer-1)){
           if(wflag[i][j][k  ]>-1000)
             wflag[i][j][k  ]=layer*int(copysign(1.0,wflag[i  ][j][k]));
        } else if(int(abs(wflag[i][j][k+1]))==ifmax &&
                  int(abs(wflag[i][j][k  ]))==(layer-1)){
           if(wflag[i][j][k+1]>-1000)
             wflag[i][j][k+1]=layer*int(copysign(1.0,wflag[i][j][k+1]));
        }
      }
    }
    insert_bc_flag(wflag, true);
    wflag.exchange();
  }
  //boil::plot->plot(clr, wflag, "clr-wflag", time->current_step());

  return;
}
