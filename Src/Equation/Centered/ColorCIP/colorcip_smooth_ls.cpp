#include "colorcip.h"
#include <iomanip>

/******************************************************************************/
void ColorCIP::smooth_ls(const Scalar & sca, Scalar & scb, const int itnum) {
/*-----------------------------------------------------------------------------+
|  convert phi of which shape is hyperbolic tangent to level set function.     |
|  reference:g.,russo and p.,smereka,"a remark on computing distance Functions"|
|            ,j.comp.physics,vol.163,2000,pp.51-67                             |
|   input    :sca, itnum                                                       |
|   output   :scb                                                              |
|   temporary:stmp                                                             |
|   flagging :iflag                                                            |
+-----------------------------------------------------------------------------*/

  const real nlayer= 6; /* Nr. layer to be calculated around the free-surface */
  const real phimin=-dxmin * nlayer;
  const real phimax= dxmin * nlayer;

  /* convert color-function [0:1] to phi [phimin:phimax] */
  for_aijk(i,j,k)
    scb[i][j][k]=(phimax-phimin)*sca[i][j][k]+phimin;

  /* store initial condition of phi */
  for_aijk(i,j,k)
    stmp[i][j][k]=scb[i][j][k];

  /*-----------+
  |  set flag  |
  +-----------*/
  for_aijk(i,j,k) {
    if(scb[i][j][k]<0.0){
      iflag[i][j][k]=-10;
    } else {
      iflag[i][j][k]= 10;
    }
  }

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=0; i<ni()-1; i++){
    for_jk(j,k){
      if(scb[i][j][k]*scb[i+1][j][k]<=0.0){
         iflag[i  ][j][k]=0;
         iflag[i+1][j][k]=0;
      }
    }
  }
  /* j-direction */
  for(int j=0; j<nj()-1; j++){
    for_ik(i,k){
      if(scb[i][j][k]*scb[i][j+1][k]<=0.0){
        iflag[i][j  ][k]=0;
        iflag[i][j+1][k]=0;
      }
    }
  }
  /* k-direction */
  for(int k=0; k<nk()-1; k++){
    for_ij(i,j){
      if(scb[i][j][k]*scb[i][j][k+1]<=0.0){
        iflag[i][j][k  ]=0;
        iflag[i][j][k+1]=0;
      }
    }
  }

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
        if(abs(iflag[i][j][k])==10&&abs(iflag[i+1][j][k])==(layer-1)){
           iflag[i  ][j][k]=layer*int(copysign(1.0,iflag[i  ][j][k]));
        } else if(abs(iflag[i+1][j][k])==10&&abs(iflag[i][j][k])==(layer-1)){
           iflag[i+1][j][k]=layer*int(copysign(1.0,iflag[i+1][j][k]));
        }
      }
    }
    /* j-direction */
    for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
        if(abs(iflag[i][j][k])==10&&abs(iflag[i][j+1][k])==(layer-1)){
           iflag[i  ][j][k]=layer*int(copysign(1.0,iflag[i  ][j][k]));
        } else if(abs(iflag[i][j+1][k])==10&&abs(iflag[i][j][k])==(layer-1)){
           iflag[i][j+1][k]=layer*int(copysign(1.0,iflag[i][j+1][k]));
        }
      }
    }
    /* k-direction */
    for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
        if(abs(iflag[i][j][k])==10&&abs(iflag[i][j][k+1])==(layer-1)){
           iflag[i  ][j][k]=layer*int(copysign(1.0,iflag[i  ][j][k]));
        } else if(abs(iflag[i][j][k+1])==10&&abs(iflag[i][j][k])==(layer-1)){
           iflag[i][j][k+1]=layer*int(copysign(1.0,iflag[i][j][k+1]));
        }
      }
    }
  }

  for(int it=1; it<=itnum; it++){
#if 1
    /* symmetric loop */
    int ist,ied,iinc,jst,jed,jinc,kst,ked,kinc;
    if((it%2)==1){
      ist=si(); ied=ei(); iinc=1;
    } else {
      ist=-ei(); ied=-si(); iinc=1;
    }
    if((it%4)==1||(it%4)==2){
      jst=sj(); jed=ej(); jinc=1;
    } else {
      jst=-ej(); jed=-sj(); jinc=1;
    }
    if(1<=(it%8)&&(it%8)<=4){
      kst=sk(); ked=ek(); kinc=1;
    } else {
      kst=-ek(); ked=-sk(); kinc=1;
    }
    for(int ic=ist; ic<=ied; ic++){
      int i=abs(ic);
    for(int jc=jst; jc<=jed; jc++){
      int j=abs(jc);
    for(int kc=kst; kc<=ked; kc++){
      int k=abs(kc);
#else
    /* normal loop */
    for_ijk(i,j,k){
#endif
      if(abs(iflag[i][j][k])<=nlayer){
        real dx = std::min(phi.dzc(k),std::min(phi.dxc(i),phi.dyc(j))); // crude
        real dtau = 0.5*dx;
        real asgn=copysign(1.0, stmp[i][j][k]);
        if(iflag[i][j][k]==0){
          /* next to free-surfacae (NFCell) */
          real dpdx1 = sqrt(
             pow((stmp[i+1][j][k]-stmp[i-1][j][k])/(phi.dxw(i)+phi.dxe(i)),2.0)
            +pow((stmp[i][j+1][k]-stmp[i][j-1][k])/(phi.dys(j)+phi.dyn(j)),2.0)
            +pow((stmp[i][j][k+1]-stmp[i][j][k-1])/(phi.dzb(k)+phi.dzt(k)),2.0));
          real dpdx2 = sqrt(
             pow((stmp[i+1][j][k]-stmp[i][j][k])/phi.dxe(i),2.0)
            +pow((stmp[i][j+1][k]-stmp[i][j][k])/phi.dyn(j),2.0)
            +pow((stmp[i][j][k+1]-stmp[i][j][k])/phi.dzt(k),2.0));
          real dpdx3 = sqrt(
             pow((stmp[i][j][k]-stmp[i-1][j][k])/phi.dxw(i),2.0)
            +pow((stmp[i][j][k]-stmp[i][j-1][k])/phi.dys(j),2.0)
            +pow((stmp[i][j][k]-stmp[i][j][k-1])/phi.dzb(k),2.0));

          real dpdx=std::max(std::max(std::max(dpdx1,dpdx2),dpdx3),1.0e-12);
          real d =  stmp[i][j][k] / dpdx;
          scb[i][j][k] -= dtau /dx * ( asgn * fabs(scb[i][j][k]) - d);
        } else {
          /* cells except NFCell */
          real a=(scb[i][j][k]-scb[i-1][j][k])/phi.dxw(i);
          real b=(scb[i+1][j][k]-scb[i][j][k])/phi.dxe(i);
          real c=(scb[i][j][k]-scb[i][j-1][k])/phi.dys(j);
          real d=(scb[i][j+1][k]-scb[i][j][k])/phi.dyn(j);
          real e=(scb[i][j][k]-scb[i][j][k-1])/phi.dzb(k);
          real f=(scb[i][j][k+1]-scb[i][j][k])/phi.dzt(k);
          real g;
          if(scb[i][j][k]>=0){
            real ap = std::max(a,0.0);
            real bm = std::min(0.0,b);
            real cp = std::max(c,0.0);
            real dm = std::min(0.0,d);
            real ep = std::max(e,0.0);
            real fm = std::min(0.0,f);
            g = sqrt(std::max(ap*ap,bm*bm)
                    +std::max(cp*cp,dm*dm)
                    +std::max(ep*ep,fm*fm))-1.0;
          } else {
            real am = std::min(0.0,a);
            real bp = std::max(b,0.0);
            real cm = std::min(0.0,c);
            real dp = std::max(d,0.0);
            real em = std::min(0.0,e);
            real fp = std::max(f,0.0);
            g = sqrt(std::max(am*am,bp*bp)
                    +std::max(cm*cm,dp*dp)
                    +std::max(em*em,fp*fp))-1.0;
          }
          scb[i][j][k] -= dtau * asgn * g; 
        }
      }
    }}}
    insert_bc_ls(scb);
    scb.exchange_all();
  }

  for_aijk(i,j,k){
    if(iflag[i][j][k]==10){
      scb[i][j][k]= 1e3;
    } else if(iflag[i][j][k]==-10){
      scb[i][j][k]=-1e3;
    }
  }
  insert_bc_ls(scb);
  scb.exchange_all();

  return;
}
