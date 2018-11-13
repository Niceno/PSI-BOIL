#include "phasechange.h"
#include <iomanip>

/******************************************************************************/
void PhaseChange::distfunc(const Scalar & sca, const int itnum) {
/***************************************************************************//**
*  /brief Convert color function to distance function.
*         Reference:G.,Russo and P.,Smereka,"A remark on computing distance
*         functions",J.Comp.Phys.,Vol.163,2000,pp.51-67
*           input    :sca, itnum
*           output   :dist
*           temporary:stmp
*           flagging :dflag
*******************************************************************************/

  //nlayer= 12; /* Nr. layer to be calculated around the free-surface */
  nlayer= 16; /* Nr. layer to be calculated around the free-surface */
  const real phimin1=-dxmin * nlayer;
  const real phimax1= dxmin * nlayer;
  const int ifmax=nlayer+4;
  const int ifmin=-nlayer-4;

  /*-------------+
  |  initialize  |
  +-------------*/
  for_aijk(i,j,k){
    delta[i][j][k]=0.0;
  }

  /*-----------+
  |  set flag  |
  +-----------*/
  for_aijk(i,j,k) {
    if(sca[i][j][k]<phisurf){         // fluid0
      dflag[i][j][k]=ifmin;
    } else {                          // fluid1
      dflag[i][j][k]=ifmax;
    }
  }

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=0; i<ni()-1; i++){
    for_jk(j,k){
      if((sca[i][j][k]-phisurf)*(sca[i+1][j][k]-phisurf)<=0.0){
         dflag[i  ][j][k]=0;
         dflag[i+1][j][k]=0;
      }
    }
  }
  /* j-direction */
  for(int j=0; j<nj()-1; j++){
    for_ik(i,k){
      if((sca[i][j][k]-phisurf)*(sca[i][j+1][k]-phisurf)<=0.0){
        dflag[i][j  ][k]=0;
        dflag[i][j+1][k]=0;
      }
    }
  }
  /* k-direction */
  for(int k=0; k<nk()-1; k++){
    for_ij(i,j){
      if((sca[i][j][k]-phisurf)*(sca[i][j][k+1]-phisurf)<=0.0){
        dflag[i][j][k  ]=0;
        dflag[i][j][k+1]=0;
      }
    }
  }
  insert_bc_flag(dflag);
  dflag.exchange();

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
        if(int(abs(dflag[i  ][j][k]))==ifmax &&
           int(abs(dflag[i+1][j][k]))==(layer-1)){
          dflag[i  ][j][k]=layer*int(copysign(1.0,dflag[i  ][j][k]));
        } else if(int(abs(dflag[i+1][j][k]))==ifmax &&
                  int(abs(dflag[i  ][j][k]))==(layer-1)){
          dflag[i+1][j][k]=layer*int(copysign(1.0,dflag[i+1][j][k]));
        }
      }
    }
    /* j-direction */
    for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
        if(int(abs(dflag[i][j  ][k]))==ifmax &&
           int(abs(dflag[i][j+1][k]))==(layer-1)){
          dflag[i  ][j][k]=layer*int(copysign(1.0,dflag[i  ][j][k]));
        } else if(int(abs(dflag[i][j+1][k]))==ifmax &&
                  int(abs(dflag[i][j  ][k]))==(layer-1)){
          dflag[i][j+1][k]=layer*int(copysign(1.0,dflag[i][j+1][k]));
        }
      }
    }
    /* k-direction */
    for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
        if(int(abs(dflag[i][j][k  ]))==ifmax &&
           int(abs(dflag[i][j][k+1]))==(layer-1)){
           dflag[i  ][j][k]=layer*int(copysign(1.0,dflag[i  ][j][k]));
        } else if(int(abs(dflag[i][j][k+1]))==ifmax &&
                  int(abs(dflag[i][j][k  ]))==(layer-1)){
          dflag[i][j][k+1]=layer*int(copysign(1.0,dflag[i][j][k+1]));
        }
      }
    }
    insert_bc_flag(dflag);
    dflag.exchange();
  }

  /*-----------------------------------------+
  |  initial condition of distance function  |
  +-----------------------------------------*/
  /* convert color-function [0:1] to phi [phimin1:phimax1] */
  for_aijk(i,j,k)
    dist[i][j][k]=(phimax1-phimin1)*sca[i][j][k]+phimin1;

  /* store initial distance function next to free-surface  */
  for_aijk(i,j,k){
    stmp[i][j][k]=(phimax1-phimin1)*sca[i][j][k]+phimin1;
  }

  if(dom->ibody().nccells() > 0) {
    for_aijk(i,j,k) {
      if(dom->ibody().off(i,j,k)){
        dflag[i][j][k]=-1001;
      }
    }
  }

  for(int it=1; it<=itnum; it++){
    for_ijk(i,j,k){
      if(abs(dflag[i][j][k])<=nlayer){  // only near region to interface
        real dx = boil::minr(phi.dxc(i),phi.dyc(j),phi.dzc(k)); // crude
        real dtau = 0.5*dx;
        real asgn=copysign(1.0, stmp[i][j][k]);
        if(dflag[i][j][k]==0){
          /* cells next to free surface (NFCell) */
          /* x-direction */
          real dpdx=0.0;
          if(dflag[i-1][j][k]==0 && dflag[i+1][j][k]==0) {
            dpdx=(stmp[i+1][j][k]-stmp[i-1][j][k])/(phi.dxw(i)+phi.dxe(i));
          } else if(dflag[i-1][j][k]==0){
            dpdx=(stmp[i][j][k]-stmp[i-1][j][k])/phi.dxw(i);
          } else if (dflag[i+1][j][k]==0) {
            dpdx=(stmp[i+1][j][k]-stmp[i][j][k])/phi.dxe(i);
          }
          /* y-direction */
          real dpdy=0.0;
          if(dflag[i][j-1][k]==0 && dflag[i][j+1][k]==0) {
            dpdy=(stmp[i][j+1][k]-stmp[i][j-1][k])/(phi.dys(j)+phi.dyn(j));
          } else if(dflag[i][j-1][k]==0){
            dpdy=(stmp[i][j][k]-stmp[i][j-1][k])/phi.dys(j);
          } else if (dflag[i][j+1][k]==0) {
            dpdy=(stmp[i][j+1][k]-stmp[i][j][k])/phi.dyn(j);
          }
          /* z-direction */
          real dpdz=0.0;
          if(dflag[i][j][k-1]==0 && dflag[i][j][k+1]==0) {
            dpdz=(stmp[i][j][k+1]-stmp[i][j][k-1])/(phi.dzb(k)+phi.dzt(k));
          } else if(dflag[i][j][k-1]==0){
            dpdz=(stmp[i][j][k]-stmp[i][j][k-1])/phi.dzb(k);
          } else if (dflag[i][j][k+1]==0) {
            dpdz=(stmp[i][j][k+1]-stmp[i][j][k])/phi.dzt(k);
          }

          real dpdxyz=sqrt(dpdx*dpdx+dpdy*dpdy+dpdz*dpdz);
          dpdxyz = std::max(dpdxyz,1.0e-12);
          real d = stmp[i][j][k] / dpdxyz;
          delta[i][j][k] = dtau /dx * ( asgn * fabs(dist[i][j][k]) - d);
        } else {
          /* cells except NFCell */
          real a=(dist[i][j][k]-dist[i-1][j][k])/phi.dxw(i);
          real b=(dist[i+1][j][k]-dist[i][j][k])/phi.dxe(i);
          real c=(dist[i][j][k]-dist[i][j-1][k])/phi.dys(j);
          real d=(dist[i][j+1][k]-dist[i][j][k])/phi.dyn(j);
          real e=(dist[i][j][k]-dist[i][j][k-1])/phi.dzb(k);
          real f=(dist[i][j][k+1]-dist[i][j][k])/phi.dzt(k);
          real g;

	  if(dflag[i-1][j][k]<-1000)a=0;
	  if(dflag[i+1][j][k]<-1000)b=0;
	  if(dflag[i][j-1][k]<-1000)c=0;
	  if(dflag[i][j+1][k]<-1000)d=0;
	  if(dflag[i][j][k-1]<-1000)e=0;
	  if(dflag[i][j][k+1]<-1000)f=0;

	  if(dflag[i][j][k]<-1000){
            g=0;
	  } else if(dist[i][j][k]>=0){
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
          delta[i][j][k] = dtau * asgn * g;
        }
      }
    }

    /* update Jakobi */
    for_ijk(i,j,k){
      if(abs(dflag[i][j][k])<=nlayer){
        dist[i][j][k] -= delta[i][j][k];
      }
    }

    insert_bc_dist(dist);
    dist.exchange_all();
  }

  for_aijk(i,j,k){
    if(dflag[i][j][k]==ifmax){
      dist[i][j][k]= real(ifmax)*dxmin;
    } else if(dflag[i][j][k]==-ifmax){
      dist[i][j][k]= real(ifmin)*dxmin;
    }
  }
  insert_bc_dist(dist);
  dist.exchange_all();

#if 0
  boil::plot->plot(clr,dist,dflag,
    "distfunc_clr-dist-dflag", time->current_step());
  exit(0);
#endif

  return;
}

