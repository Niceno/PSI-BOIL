#include "cipcsl2.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void CIPCSL2::ext_sca(Scalar & sca) {
/***************************************************************************//**
*  \brief extrapolate sca 
*         output : sca
*         temporary: iflag,fn
*******************************************************************************/

#ifdef DEBUG
  boil::oout<<"cipcsl2_ext_fs:start\n";
#endif
  for_avijk(iflag,i,j,k) {
    iflag[i][j][k]=0;
  }
  /* i-direction */
  for(int i=sca.si()-1; i<=sca.ei(); i++){
    for_vjk(sca,j,k){
       if((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0){
          iflag[i  ][j][k]=1;
          iflag[i+1][j][k]=1;
       }
    }
  }
  /* j-direction */
  for(int j=sca.sj()-1; j<=sca.ej(); j++){
    for_vik(sca,i,k){
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
          iflag[i][j  ][k]=1;
          iflag[i][j+1][k]=1;
       }
    }
  }
  /* k-direction */
  for(int k=sca.sk()-1; k<=sca.ek(); k++){
    for_vij(sca,i,j){
       if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
          iflag[i][j][k  ]=1;
          iflag[i][j][k+1]=1;
       }
    }
  }
  iflag.exchange_all();

  /*-------------+
  |  advect alp  |
  +-------------*/

  //boil::plot->plot(clr,dist,iflag,stmp, 
  //                "clr-dist-iflag-stmp", time->current_step());
  //exit(0);

  int mmax=4;

  for(int mstep=1; mstep<=mmax; mstep++){
    for(int it=1; it<=4; it++){
      int ist,ied,iinc;
      int jst,jed,jinc;
      int kst,ked,kinc;
      if(it%2==0){ist=si();ied=ei();iinc=1;}else{ist=ei();ied=si();iinc=-1;}
      if(it%2==0){jst=sj();jed=ej();jinc=1;}else{jst=ej();jed=sj();jinc=-1;}
      if(it%2==0){kst=sk();ked=ek();kinc=1;}else{kst=ek();ked=sk();kinc=-1;}
      for(int i=ist; i<=ied; i+=iinc){
      for(int j=jst; j<=jed; j+=jinc){
      for(int k=kst; k<=ked; k+=kinc){
        if(int(iflag[i][j][k])==0) {
          real dalpsum=0.0;
          int isum=0;
          if(int(iflag[i-1][j][k])==1){
            dalpsum+=sca[i-1][j][k];
            isum++;
          }
          if(int(iflag[i+1][j][k])==1){
            dalpsum+=sca[i+1][j][k];
            isum++;
          }
          if(int(iflag[i][j-1][k])==1){
            dalpsum+=sca[i][j-1][k];
            isum++;
          }
          if(int(iflag[i][j+1][k])==1){
            dalpsum+=sca[i][j+1][k];
            isum++;
          }
          if(int(iflag[i][j][k-1])==1){
            dalpsum+=sca[i][j][k-1];
            isum++;
          }
          if(int(iflag[i][j][k+1])==1){
            dalpsum+=sca[i][j][k+1];
            isum++;
          }
          fn[i][j][k]=dalpsum/real(isum+1.0e-12);
        }
      }
    }}}

    /* update sca */
    for_ijk(i,j,k) {
      if(int(iflag[i][j][k])==0){
        sca[i][j][k] = fn[i][j][k];
      }
    }
    sca.exchange_all();
  }

#ifdef DEBUG
  boil::oout<<"cipcsl2_set_alp:end\n";
#endif

  return;
}
