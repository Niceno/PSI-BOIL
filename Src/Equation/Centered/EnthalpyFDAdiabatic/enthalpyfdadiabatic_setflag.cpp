#include "enthalpyfdadiabatic.h"

/***************************************************************************//**
*  \brief Set flag for convection term
*  iflag = -1 Finite Volume Method. color function <  clrsurf
*  iflag = 1  Finite Volume Method. color function >= clrsurf
*  iflag = 0  Finite Difference Method
*  iflag = -1000 solid cell (Finite Volume Method)
*******************************************************************************/
void EnthalpyFDAdiabatic::setflag() {

  for_aijk(i,j,k) {
    iflag[i][j][k]=-1;
  }
  for_aijk(i,j,k) {
    if(clrold[i][j][k]>=clrsurf)
      iflag[i][j][k]=1;
  }
  /* i-direction */
  for(int i=si()-1; i<=ei(); i++){
    for_jk(j,k){
       if((clrold[i][j][k]-clrsurf)*(clrold[i+1][j][k]-clrsurf)<=0.0){
          iflag[i  ][j][k]=0;
          iflag[i+1][j][k]=0;
       }
    }
  }
  /* j-direction */
  for(int j=sj()-1; j<=ej(); j++){
    for_ik(i,k){
      if((clrold[i][j][k]-clrsurf)*(clrold[i][j+1][k]-clrsurf)<=0.0){
          iflag[i][j  ][k]=0;
          iflag[i][j+1][k]=0;
       }
    }
  }
  /* k-direction */
  for(int k=sk()-1; k<=ek(); k++){
    for_ij(i,j){
       if((clrold[i][j][k]-clrsurf)*(clrold[i][j][k+1]-clrsurf)<=0.0){
          iflag[i][j][k  ]=0;
          iflag[i][j][k+1]=0;
       }
    }
  }

  for_ijk(i,j,k) {
    if (dom->ibody().off(i,j,k)) 
      iflag[i][j][k]=-1000;
  }

  iflag.exchange();
}
