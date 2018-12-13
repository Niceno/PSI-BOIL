#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::cal_gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate gradient of temperature
*         txl,tyl,tzl: gradient of liquid tempereture in x,y,z-direction
*         txv,tyv,tzv: gradient of vapor  tempereture in x,y,z-direction
*******************************************************************************/

  /* this is for increased precision near cell centre */
  /* achieved by switching to upwind/downwind differencing */
  /* it is also possible to use purely upw/dwnw through the upwind_flag */
  /* this is however not recommended <- lower precision */
  prepare_gradt8();

  for_vijk(tpr,i,j,k){

    /* normal cell, 2nd order */
    real dtdx,dtdy,dtdz;
    tpr.grad(i,j,k,&dtdx,&dtdy,&dtdz);
    if(clr[i][j][k]>=clrsurf){
      txl[i][j][k]=dtdx;
      tyl[i][j][k]=dtdy;
      tzl[i][j][k]=dtdz;
    } else {
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

    /* the clrsurf check should work as long as update_at_walls
     * is called in VOF! */

    /* west */
    if((clrw-clrsurf)*(clrc-clrsurf)<=0.0){
      txp = gradtx(-1,i,j,k);
      if((clrc-clrsurf)<0.0){
        txv[i][j][k] = txp; 
      } else {
        txl[i][j][k] = txp;
      }
      ii=1;
    }
    /* east */
    if((clre-clrsurf)*(clrc-clrsurf)<=0.0){
      txm = gradtx(+1,i,j,k);
      if((clrc-clrsurf)<0.0){
        txv[i][j][k]=txm;
      } else {
        txl[i][j][k]=txm;
      }
      ii+=1;
    }
    /* south */
    if((clrs-clrsurf)*(clrc-clrsurf)<=0.0){
      typ = gradty(-1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tyv[i][j][k]=typ;
      } else {
        tyl[i][j][k]=typ;
      }
      jj=1;
    }
    /* north */
    if((clrn-clrsurf)*(clrc-clrsurf)<=0.0){
      tym = gradty(+1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tyv[i][j][k]=tym;
      } else {
        tyl[i][j][k]=tym;
      }
      jj+=1;
    }
    /* bottom */
    if((clrb-clrsurf)*(clrc-clrsurf)<=0.0){
      tzp = gradtz(-1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tzv[i][j][k]=tzp;
      } else {
        tzl[i][j][k]=tzp;
      }
      kk=1;
    }
    /* top */
    if((clrt-clrsurf)*(clrc-clrsurf)<=0.0){
      tzm = gradtz(+1,i,j,k);
      if((clrc-clrsurf)<0.0){
        tzv[i][j][k]=tzm;
      } else {
        tzl[i][j][k]=tzm;
      }
      kk+=1;
    }

#if 1
    if(ii==2)txv[i][j][k]=txl[i][j][k]=0.0;
    if(jj==2)tyv[i][j][k]=tyl[i][j][k]=0.0;
    if(kk==2)tzv[i][j][k]=tzl[i][j][k]=0.0;
#endif
  }

#ifdef IB
  gradt_ib(diff_eddy);
#endif
  
  /* correct at walls */
  insert_bc_gradt(diff_eddy);

  txl.exchange_all();
  tyl.exchange_all();
  tzl.exchange_all();
  txv.exchange_all();
  tyv.exchange_all();
  tzv.exchange_all();

  return;
}
