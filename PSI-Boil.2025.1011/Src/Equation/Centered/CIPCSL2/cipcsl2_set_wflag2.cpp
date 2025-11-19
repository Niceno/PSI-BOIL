#include "cipcsl2.h"
#include <iomanip>
using namespace std;

/******************************************************************************/
void CIPCSL2::set_wflag2() {
/***************************************************************************//**
*  /brief flagging
*           output   :wflag
*
*   wflag: -3 | -2 | -1 |  0 |  0 |  1 |  2 |  3 | ..   
*          ---+----+----+----+----+----+----+----+---
*   vapor                   |                        liquid
*                      free-surface
*
*   solid region                   : wflag = -1001
*   cell except wall adjacent cells: wflag = -1001
*******************************************************************************/
#ifdef DEBUG
  cout<<"set_wflag2: begin\n";
#endif
  const int ifmax =  nlayer+4;
  const int ifmin = -nlayer-4;

  /*-------------+
  |  reset flag  |
  +-------------*/
  wflag = -1001;

  /*----------------------------+
  |  wall adjacent cell = 1001  |
  |  otherwise = -1001          |
  +----------------------------*/
  for( int b=0; b<phi.bc().count(); b++ ) {
    if(phi.bc().type_decomp(b)) continue;
    if( phi.bc().type(b) == BndType::wall() ) {
      Dir d      = phi.bc().direction(b);
      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        int iof=0, jof=0, kof=0;
        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        int abs_iof = abs(iof);
        int abs_jof = abs(jof);
        int abs_kof = abs(kof);

        for_vijk( phi.bc().at(b), i,j,k ){
          if (i<si()-abs_iof) continue;
          if (i>ei()+abs_iof) continue;
          if (j<sj()-abs_jof) continue;
          if (j>ej()+abs_jof) continue;
          if (k<sk()-abs_kof) continue;
          if (k>ek()+abs_kof) continue;
          wflag[i+iof][j+jof][k+kof]=1001;
        }
      }
    }
  }

  /*-----------------------------------+
  |  set flag for wall adjacent cells  |
  +-----------------------------------*/
  for_ijk(i,j,k) {
    if (wflag[i][j][k]>1000) {
      if(clr[i][j][k]<phisurf){
        wflag[i][j][k]=ifmin;
      } else {
        wflag[i][j][k]=ifmax;
      }
    }
  }
  wflag.exchange();

#ifdef IB
  /* wall adjacent cells */
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    if(clr[i][j][k]<phisurf){
      wflag[i][j][k]=ifmin;
    } else {
      wflag[i][j][k]=ifmax;
    }
  }
#endif
  wflag.exchange();

#if 0
  cout<<"set_wflag2: "<<wflag[4][17][1]<<" "<<clr[4][17][1]<<" "<<clr[3][17][1]<<"\n";
#endif

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=si()-1; i<ei()+1; i++){
    for_jk(j,k){
      if ((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0) {
        if(wflag[i  ][j][k]>-1000 && wflag[i+1][j][k]>-1000) {
          wflag[i  ][j][k]=0;
          wflag[i+1][j][k]=0;
        }
      }
    }
  }
  /* j-direction */
  for(int j=sj()-1; j<ej()+1; j++){
    for_ik(i,k){
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
        if(wflag[i][j  ][k]>-1000 && wflag[i][j+1][k]>-1000) {
          wflag[i][j  ][k]=0;
          wflag[i][j+1][k]=0;
        }
      }
    }
  }
  /* k-direction */
  for(int k=sk()-1; k<ek()+1; k++){
    for_ij(i,j){
      if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
        if(wflag[i][j][k  ]>-1000 && wflag[i][j][k+1]>-1000) {
           wflag[i][j][k  ]=0;
           wflag[i][j][k+1]=0;
        }
      }
    }
  }
  wflag.exchange();
#if 0
  cout<<"set_wflag2: "<<wflag[4][17][1]<<" "<<clr[4][17][1]<<" "<<clr[3][17][1]<<" "<<clr[4][17][0]<<" "<<clr[4][17][2]<<"\n";
#endif


  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=si()-1; i<ei()+1; i++){
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
    for(int j=sj()-1; j<ej()+1; j++){
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
    for(int k=sk()-1; k<ek()+1; k++){
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
    wflag.exchange();
  }

#if 0
  boil::plot->plot(clr, wflag, "clr-wflag", time->current_step());
  exit(0);
#endif

  return;
}
