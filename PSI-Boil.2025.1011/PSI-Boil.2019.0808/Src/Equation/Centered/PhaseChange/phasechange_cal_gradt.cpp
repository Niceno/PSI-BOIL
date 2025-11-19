#include "phasechange.h"
using namespace std;

real temperature_node(real len_s, real lam_s, real tmp_s
                    , real len_f, real lam_f, real tmp_f);
real grad3(real ww, real dm, real dp, real tm, real tc, real tp, real epsl);
/******************************************************************************/
void PhaseChange::cal_gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate gradient of temperature
*         txl,tyl,tzl: gradient of liquid tempereture in x,y,z-direction
*         txv,tyv,tzv: gradient of vapor  tempereture in x,y,z-direction
*******************************************************************************/

  for_vijk(tpr,i,j,k){

    /* normal cell, 2nd order */
    real dtdx,dtdy,dtdz;
    tpr.grad(i,j,k,&dtdx,&dtdy,&dtdz);
    if(clr[i][j][k]>=phisurf){
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

    /* west */
    if((clrw-phisurf)*(clrc-phisurf)<=0.0){
      gradtx5( i-1, j, k, &txm, &txp, 1);
      if((clrc-phisurf)<0.0){
        txv[i][j][k]=txp;
      } else {
        txl[i][j][k]=txp;
      }
      ii=1;
    }
    /* east */
    if((clre-phisurf)*(clrc-phisurf)<=0.0){
      gradtx5( i, j, k, &txm, &txp, -1);
      if((clrc-phisurf)<0.0){
        txv[i][j][k]=txm;
      } else {
        txl[i][j][k]=txm;
      }
      ii+=1;
    }
    /* south */
    if((clrs-phisurf)*(clrc-phisurf)<=0.0){
      gradty5( i, j-1, k, &tym, &typ, 1);
      if((clrc-phisurf)<0.0){
        tyv[i][j][k]=typ;
      } else {
        tyl[i][j][k]=typ;
      }
      jj=1;
    }
    /* north */
    if((clrn-phisurf)*(clrc-phisurf)<=0.0){
      gradty5( i, j, k, &tym, &typ, -1);
      if((clrc-phisurf)<0.0){
        tyv[i][j][k]=tym;
      } else {
        tyl[i][j][k]=tym;
      }
      jj+=1;
    }
    /* bottom */
    if((clrb-phisurf)*(clrc-phisurf)<=0.0){
      gradtz5( i, j, k-1, &tzm, &tzp, 1);
      if((clrc-phisurf)<0.0){
        tzv[i][j][k]=tzp;
      } else {
        tzl[i][j][k]=tzp;
      }
      kk=1;
    }
    /* top */
    if((clrt-phisurf)*(clrc-phisurf)<=0.0){
      gradtz5( i, j, k, &tzm, &tzp, -1);
      if((clrc-phisurf)<0.0){
        tzv[i][j][k]=tzm;
      } else {
        tzl[i][j][k]=tzm;
      }
      kk+=1;
    }
#if 0
#else
    if(ii==2)txv[i][j][k]=txl[i][j][k]=0.0;
    if(jj==2)tyv[i][j][k]=tyl[i][j][k]=0.0;
    if(kk==2)tzv[i][j][k]=tzl[i][j][k]=0.0;
#endif
  }

#ifdef IB
  /*--------------------------+ 
  | touch from immersed body  |
  +--------------------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    real t_c  = tpr[i][j][k];

    /* x direction */
    Comp m = Comp::u();
    real dx_w = phi.dxw(i);
    real dx_e = phi.dxe(i);
    real t_w  = tpr[i-1][j][k];
    real t_e  = tpr[i+1][j][k];

    // west is in wall
    if(dom->ibody().off(m,i,j,k)) {
      dx_w *= dom->ibody().fdxw(cc);
      real len_s = phi.dxw(i) - 0.5*phi.dxc(i);
      real lam_s = solid()->lambda(i-1,j,k);
      real tmp_s = tpr[i-1][j][k];
      real len_f = 0.5*phi.dxc(i);
      real lam_f;
      if(clr[i][j][k]<phisurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
      real dtdx = (t_e-t_w)/(dx_w+dx_e);

      // interface in east
      if((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0){
        t_e = tsat;
        real ww = (phisurf-clr[i][j][k])/(clr[i+1][j][k]-clr[i][j][k]);
        dx_e *= ww;
	dtdx = grad3(ww, dx_w, dx_e, t_w, t_c, t_e, epsl);
      }

      // update grad
      if(clr[i][j][k]>=phisurf) {
        txl[i][j][k]=dtdx;
      } else {
        txv[i][j][k]=dtdx;
      }
    }

    // east is in wall
    if(dom->ibody().off(m,i+1,j,k)) {
      dx_e *= dom->ibody().fdxe(cc);
      real len_s = phi.dxe(i) - 0.5*phi.dxc(i);
      real lam_s = solid()->lambda(i+1,j,k);
      real tmp_s = tpr[i+1][j][k];
      real len_f = 0.5*phi.dxc(i);
      real lam_f;
      if(clr[i][j][k]<phisurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_e = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
      real dtdx = (t_e-t_w)/(dx_w+dx_e);

      // interface in west
      if((clr[i-1][j][k]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){
        t_w = tsat;
        real ww = (phisurf-clr[i][j][k])/(clr[i-1][j][k]-clr[i][j][k]);
        dx_w *= ww;
	dtdx = grad3(ww, dx_w, dx_e, t_w, t_c, t_e, epsl);
      }

      // update grad
      if(clr[i][j][k]>=phisurf) {
        txl[i][j][k]=dtdx;
      } else {
        txv[i][j][k]=dtdx;
      }
    }

    /* y direction */
    m = Comp::v();
    real dy_s = phi.dys(j);
    real dy_n = phi.dyn(j);
    real t_s  = tpr[i][j-1][k];
    real t_n  = tpr[i][j+1][k];

    // south is in wall
    if(dom->ibody().off(m,i,j,k)) {
      dy_s *= dom->ibody().fdys(cc);
      real len_s = phi.dys(j) - 0.5*phi.dyc(j);
      real lam_s = solid()->lambda(i,j-1,k);
      real tmp_s = tpr[i][j-1][k];
      real len_f = 0.5*phi.dyc(j);
      real lam_f;
      if(clr[i][j][k]<phisurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_s = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
      real dtdy = (t_n-t_s)/(dy_s+dy_n);

      // interface in north
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
        t_n = tsat;
        real ww = (phisurf-clr[i][j][k])/(clr[i][j+1][k]-clr[i][j][k]);
        dy_n *= ww;
	dtdy = grad3(ww, dy_s, dy_n, t_s, t_c, t_n, epsl);
      }

      // update grad
      if(clr[i][j][k]>=phisurf) {
        tyl[i][j][k]=dtdy;
      } else {
        tyv[i][j][k]=dtdy;
      }
    }

    // north is in wall
    if(dom->ibody().off(m,i,j+1,k)) {
      dy_n *= dom->ibody().fdyn(cc);
      real len_s = phi.dyn(j) - 0.5*phi.dyc(j);
      real lam_s = solid()->lambda(i,j+1,k);
      real tmp_s = tpr[i][j+1][k];
      real len_f = 0.5*phi.dyc(j);
      real lam_f;
      if(clr[i][j][k]<phisurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_n = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
      real dtdy = (t_n-t_s)/(dy_s+dy_n);

      // interface in south
      if((clr[i][j-1][k]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){
        t_s = tsat;
        real ww = (phisurf-clr[i][j][k])/(clr[i][j-1][k]-clr[i][j][k]);
        dy_s *= ww;
	dtdy = grad3(ww, dy_s, dy_n, t_s, t_c, t_n, epsl);
      }

      // update grad
      if(clr[i][j][k]>=phisurf) {
        tyl[i][j][k]=dtdy;
      } else {
        tyv[i][j][k]=dtdy;
      }
    }

    // z direction
    m = Comp::w();
    real dz_b = phi.dzb(k);
    real dz_t = phi.dzt(k);
    real t_b  = tpr[i][j][k-1];
    real t_t  = tpr[i][j][k+1];

    // bottom is in wall
    if(dom->ibody().off(m,i,j,k)) {
      dz_b *= dom->ibody().fdzb(cc);
#if 1
      real len_s = phi.dzb(k) - 0.5*phi.dzc(k);
      real lam_s = solid()->lambda(i,j,k-1);
      real tmp_s = tpr[i][j][k-1];
      real len_f = 0.5*phi.dzc(k);
      real lam_f;
      if(clr[i][j][k]<phisurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_b = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_b = 109.444;  // special setting for verification!!!
#endif
      real dtdz = (t_t-t_b)/(dz_b+dz_t);

      // interface in top
      if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
        t_t = tsat;
	real ww = (phisurf-clr[i][j][k])/(clr[i][j][k+1]-clr[i][j][k]);
       	dz_t *= ww;
	dtdz = grad3(ww, dz_b, dz_t, t_b, t_c, t_t, epsl);
        if (ww>epsl) {
          real a = dz_b;
          real b = dz_t;
          dtdz = b*b*(t_c - t_b) +a*a*(t_t - t_c);
          dtdz /= (a*b*(a+b));
        } else {
          dtdz = (t_t - t_b)/(dz_b + dz_t);
        }
      }

      // update grad
      if(clr[i][j][k]>=phisurf) {
        tzl[i][j][k]=dtdz;
      } else {
        tzv[i][j][k]=dtdz;
      }
    }

    // top is in wall
    if(dom->ibody().off(m,i,j,k+1)) {
      dz_t *= dom->ibody().fdzt(cc);
      real len_s = phi.dzt(k) - 0.5*phi.dzc(k);
      real lam_s = solid()->lambda(i,j,k+1);
      real tmp_s = tpr[i][j][k+1];
      real len_f = 0.5*phi.dzc(k);
      real lam_f;
      if(clr[i][j][k]<phisurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_t = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
      real dtdz = (t_t-t_b)/(dz_b+dz_t);

      // interface in bottom
      if((clr[i][j][k-1]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){
        t_b = tsat;
        real ww = (phisurf-clr[i][j][k])/(clr[i][j][k-1]-clr[i][j][k]);
        dz_b *= ww;
	dtdz = grad3(ww, dz_b, dz_t, t_b, t_c, t_t, epsl);
      }

      // update grad
      if(clr[i][j][k]>=phisurf) {
        tzl[i][j][k]=dtdz;
      } else {
        tzv[i][j][k]=dtdz;
      }
    }
  }

  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) {
        txl[i][j][k]=0.0;
        tyl[i][j][k]=0.0;
        tzl[i][j][k]=0.0;
        txv[i][j][k]=0.0;
        tyv[i][j][k]=0.0;
        tzv[i][j][k]=0.0;
      }
    }
  }
#endif


  txl.exchange_all();
  tyl.exchange_all();
  tzl.exchange_all();
  txv.exchange_all();
  tyv.exchange_all();
  tzv.exchange_all();

  return;
}

/******************************************************************************/
real temperature_node(real len_s, real lam_s, real tmp_s
                    , real len_f, real lam_f, real tmp_f) {
/***************************************************************************//**
*  \brief calculate temperature at node point
*             len_s         len_f
*             lam_s         lam_f
*         *-------------*------------*
*       tmp_s       tmp_node        tmp_f
*******************************************************************************/
  return (len_f*lam_s*tmp_s + len_s*lam_f*tmp_f)/(len_f*lam_s + len_s*lam_f);
}

/******************************************************************************/
real grad3(real ww, real dm, real dp, real tm, real tc, real tp, real epsl){
   real dtdx;
/***************************************************************************//**
*  \brief calculate gradient from three point
*               dm           dp
*         *-------------*------------*
*        tm            tc           tp 
*******************************************************************************/
  if (ww>epsl) {
    real a = dm;
    real b = dp;
    dtdx = b*b*(tc - tm) +a*a*(tp - tc);
    dtdx /= (a*b*(a+b));
   } else {
    dtdx = (tp - tm)/(dm + dp);
   }
  return dtdx;
}
