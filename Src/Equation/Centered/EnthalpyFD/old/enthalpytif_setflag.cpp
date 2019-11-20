#include "enthalpytif.h"

/***************************************************************************//**
*  \brief Set flag for convection term
*  iflag = -1 Finite Volume Method. color function <  clrsurf
*  iflag = 1  Finite Volume Method. color function >= clrsurf
*  iflag = 0  Finite Difference Method
*  iflag = -1000 solid cell (Finite Volume Method)
*******************************************************************************/
void EnthalpyTIF::setflag() {

  for_aijk(i,j,k) {
    iflag[i][j][k]=-1;
    if(clrold[i][j][k]>=clrsurf)
      iflag[i][j][k]=1;
  }

  /* i-direction */
  for(int i=si()-1; i<=ei(); i++){
    for_jk(j,k){
      //if((clrold[i][j][k]-clrsurf)*(clrold[i+1][j][k]-clrsurf)<=0.0
      //   &&iflag[i][j][k]+iflag[i+1][j][k]>-500) {
      if(Interface_old(+1,Comp::i(),i,j,k)) {
        iflag[i  ][j][k]=0;
        iflag[i+1][j][k]=0;
      }
    }
  }
  /* j-direction */
  for(int j=sj()-1; j<=ej(); j++){
    for_ik(i,k){
      //if((clrold[i][j][k]-clrsurf)*(clrold[i][j+1][k]-clrsurf)<=0.0
      //    &&iflag[i][j][k]+iflag[i][j+1][k]>-500) {
      if(Interface_old(+1,Comp::j(),i,j,k)) {
        iflag[i][j  ][k]=0;
        iflag[i][j+1][k]=0;
      }
    }
  }
  /* k-direction */
  for(int k=sk()-1; k<=ek(); k++){
    for_ij(i,j){
      //if((clrold[i][j][k]-clrsurf)*(clrold[i][j][k+1]-clrsurf)<=0.0
      //    &&iflag[i][j][k]+iflag[i][j][k+1]>-500) {
      if(Interface_old(+1,Comp::j(),i,j,k)) {
        iflag[i][j][k  ]=0;
        iflag[i][j][k+1]=0;
      }
    }
  }

  for_aijk(i,j,k) {
    if(dom->ibody().off(i,j,k)) 
      iflag[i][j][k]=-1000;
  }

#if 0
  /* special marching-cube based trigger for interfacial cells
   * that are only partially cut by the isosurface */
  for_ijk(i,j,k) {
    if(iflag[i][j][k] != 0 && iflag[i][j][k] != -1000) {
      real nodeval;
      for (int m=0; m<=7; m++){  
        int ii,jj,kk;
        if(m==0)     {ii=i-1; jj=j-1; kk=k-1;}
        else if(m==1){ii=i  ; jj=j-1; kk=k-1;}
        else if(m==2){ii=i  ; jj=j  ; kk=k-1;}
        else if(m==3){ii=i-1; jj=j  ; kk=k-1;}
        else if(m==4){ii=i-1; jj=j-1; kk=k  ;}
        else if(m==5){ii=i  ; jj=j-1; kk=k  ;}
        else if(m==6){ii=i  ; jj=j  ; kk=k  ;}
        else if(m==7){ii=i-1; jj=j  ; kk=k  ;}
        nodeval=(clrold[ii][jj  ][kk  ]+clrold[ii+1][jj  ][kk  ]
                +clrold[ii][jj+1][kk  ]+clrold[ii+1][jj+1][kk  ]
                +clrold[ii][jj  ][kk+1]+clrold[ii+1][jj  ][kk+1]
                +clrold[ii][jj+1][kk+1]+clrold[ii+1][jj+1][kk+1])/8.0;

        if((nodeval-clrsurf)*(clrold[i][j][k]-clrsurf)<=0.0) {
          iflag[i][j][k]=0;
          continue;
        }  
      }
    }
  }
#else
  /* trigger based on adens */
  if(adens)
    for_ijk(i,j,k) {
      if(iflag[i][j][k] != 0 && iflag[i][j][k] != -1000) {
        if(adensold[i][j][k]>boil::pico) {
          iflag[i][j][k]=0;
        } 
      }
    }
#endif

  iflag.exchange();
}	

