#include "phasechange.h"
#include <iomanip>
//#define DEBUG

/******************************************************************************/
void PhaseChange::ib_set_dflag() {
/***************************************************************************//**
*  /brief flagging for immersed boundary
*  dflag = 2: fluid far from solid
*  dflag = 1: fluid next to solid
*  dflag =-1: solid next to fluid
*  dflag =-2: solid far from fluid
*           output   :dflag
*******************************************************************************/

#ifdef DEBUG
  std::cout<<"ib_set_dflag:begin \n";
#endif

  /*-----------+
  |  set flag  |
  +-----------*/
  if(dom->ibody().nccells() > 0) {
    for_aijk(i,j,k) {
      if(dom->ibody().fV(i,j,k)<0.5){
        dflag[i][j][k]=-2;
      } else {
        dflag[i][j][k]=2;
      }
    }
  } else {
    dflag=2;
  }

#ifdef DEBUG
  std::cout<<"ib_set_dflag:set 2 or -2 \n";
#endif

  insert_bc_flag(dflag);
  dflag.exchange();

#ifdef DEBUG
  std::cout<<"ib_set_dflag:set 2 or -2 \n";
#endif

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    dflag[i][j][k]=1;
  }
  insert_bc_flag(dflag);
  dflag.exchange();

  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    dflag[i][j][k]=1.0;
    if(dflag[i-1][j][k]<0)dflag[i-1][j][k]=-1;
    if(dflag[i+1][j][k]<0)dflag[i+1][j][k]=-1;
    if(dflag[i][j-1][k]<0)dflag[i][j-1][k]=-1;
    if(dflag[i][j+1][k]<0)dflag[i][j+1][k]=-1;
    if(dflag[i][j][k-1]<0)dflag[i][j][k-1]=-1;
    if(dflag[i][j][k+1]<0)dflag[i][j][k+1]=-1;
  }
  insert_bc_flag(dflag);
  dflag.exchange();

#if 0
  boil::plot->plot(clr,dflag, "ib_set_dflag_clr-dflag", time->current_step());
  exit(0);
#endif

  return;
}
