#include "cipcsl2.h"
#include <iomanip>
//#define DEBUG

/******************************************************************************/
void CIPCSL2::ib_set_iflag() {
/***************************************************************************//**
*  /brief flagging for immersed boundary
*  iflag = 2: fluid far from solid
*  iflag = 1: fluid next to solid
*  iflag =-1: solid next to fluid
*  iflag =-2: solid far from fluid
*           output   :iflag
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"ib_set_iflag:begin \n";
#endif
  /*-----------+
  |  set flag  |
  +-----------*/
  if(dom->ibody().nccells() > 0) {
    for_aijk(i,j,k) {
      if(dom->ibody().fV(i,j,k)<0.5){
        iflag[i][j][k]=-2;
      } else {
        iflag[i][j][k]=2;
      }
    }
  } else {
    iflag=2;
  }

#ifdef DEBUG
  std::cout<<"ib_set_iflag:set 2 or -2 \n";
#endif

  insert_bc_flag(iflag, false);
  iflag.exchange();

#ifdef DEBUG
  std::cout<<"ib_set_iflag:set 2 or -2 \n";
#endif

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    iflag[i][j][k]=1;
  }
  insert_bc_flag(iflag, false);
  iflag.exchange();

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    iflag[i][j][k]=1;
    if(iflag[i-1][j][k]<0)iflag[i-1][j][k]=-1;
    if(iflag[i+1][j][k]<0)iflag[i+1][j][k]=-1;
    if(iflag[i][j-1][k]<0)iflag[i][j-1][k]=-1;
    if(iflag[i][j+1][k]<0)iflag[i][j+1][k]=-1;
    if(iflag[i][j][k-1]<0)iflag[i][j][k-1]=-1;
    if(iflag[i][j][k+1]<0)iflag[i][j][k+1]=-1;
  }
  insert_bc_flag(iflag, false);
  iflag.exchange();

  //boil::plot->plot(clr,iflag, "ib_set_iflag_clr-iflag", time->current_step());

  return;
}
