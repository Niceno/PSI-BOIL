#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::cal_gradt_fluid(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate gradient of temperature
*         txl,tyl,tzl: gradient of liquid tempereture in x,y,z-direction
*         txv,tyv,tzv: gradient of vapor  tempereture in x,y,z-direction
*******************************************************************************/

  for_vijk(tpr,i,j,k){

    if(dom->ibody().off(i,j,k)) continue;

    /* normal cell, 2nd order */
    real dtdx,dtdy,dtdz;
    tpr.grad(i,j,k,&dtdx,&dtdy,&dtdz);
    if(clr[i][j][k]>=clrsurf){
      txl[i][j][k]=dtdx;
      tyl[i][j][k]=dtdy;
      tzl[i][j][k]=dtdz;
      txv[i][j][k]=0.0;
      tyv[i][j][k]=0.0;
      tzv[i][j][k]=0.0;
    } else {
      txl[i][j][k]=0.0;
      tyl[i][j][k]=0.0;
      tzl[i][j][k]=0.0;
      txv[i][j][k]=dtdx;
      tyv[i][j][k]=dtdy;
      tzv[i][j][k]=dtdz;
    }

    /* interface cell, 2nd order */
    real clrc=clr[i][j][k];
    real clrw=clr[i-1][j][k];
    real clre=clr[i+1][j][k];
    real clrs=clr[i][j-1][k];
    real clrn=clr[i][j+1][k];
    real clrb=clr[i][j][k-1];
    real clrt=clr[i][j][k+1];
    real txm, tym, tzm, txp, typ, tzp;
    int ii,jj,kk;

    ii=jj=kk=0;

    /* west */
    //if((clrw-clrsurf)*(clrc-clrsurf)<=0.0){
    if(Interface(-1,Comp::i(),i,j,k)) {
      txp = gradtx(-1,i,j,k);
      if((clrc-clrsurf)<0.0){
        txv[i][j][k] = txp; 
        txl[i][j][k] = 0.0;
      } else {
        txl[i][j][k] = txp;
        txv[i][j][k] = 0.0;
      }
      ii=1;
    }
    /* east */
    //if((clre-clrsurf)*(clrc-clrsurf)<=0.0){
    if(Interface(+1,Comp::i(),i,j,k)) {
      txm = gradtx(+1,i,j,k);
      if((clrc-clrsurf)<0.0){
        txv[i][j][k] = txm;
        txl[i][j][k] = 0.0;
      } else {
        txl[i][j][k]=txm;
        txv[i][j][k] = 0.0;
      }
      ii+=1;
    }
    /* south */
    //if((clrs-clrsurf)*(clrc-clrsurf)<=0.0){
    if(Interface(-1,Comp::j(),i,j,k)) {
      typ = gradty(-1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tyv[i][j][k] = typ;
        tyl[i][j][k] = 0.0;
      } else {
        tyl[i][j][k]=typ;
        tyv[i][j][k] = 0.0;
      }
      jj=1;
    }
    /* north */
    //if((clrn-clrsurf)*(clrc-clrsurf)<=0.0){
    if(Interface(+1,Comp::j(),i,j,k)) {
      tym = gradty(+1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tyv[i][j][k] = tym;
        tyl[i][j][k] = 0.0;
      } else {
        tyl[i][j][k] = tym;
        tyv[i][j][k] = 0.0;
      }
      jj+=1;
    }
    /* bottom */
    //if((clrb-clrsurf)*(clrc-clrsurf)<=0.0){
    if(Interface(-1,Comp::k(),i,j,k)) {
      tzp = gradtz(-1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tzv[i][j][k] = tzp;
        tzl[i][j][k] = 0.0;
      } else {
        tzl[i][j][k] = tzp;
        tzv[i][j][k] = 0.0;
      }
      kk=1;
    }
    /* top */
    //if((clrt-clrsurf)*(clrc-clrsurf)<=0.0){
    if(Interface(+1,Comp::k(),i,j,k)) {
      tzm = gradtz(+1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tzv[i][j][k] = tzm;
        tzl[i][j][k] = 0.0;
      } else {
        tzl[i][j][k] = tzm;
        tzv[i][j][k] = 0.0;
      }
      kk+=1;
    }

#if 1
    if(ii==2)txv[i][j][k]=txl[i][j][k]=0.0;
    if(jj==2)tyv[i][j][k]=tyl[i][j][k]=0.0;
    if(kk==2)tzv[i][j][k]=tzl[i][j][k]=0.0;
#endif
  }

  return;
}
