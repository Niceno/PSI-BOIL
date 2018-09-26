#include "cipcsl2.h"
#include <iomanip>

/******************************************************************************/
void CIPCSL2::distfunc(const Scalar & sca, const int itnum) {
/***************************************************************************//**
*  /brief Convert color function to distance function.
*         Reference:G.,Russo and P.,Smereka,"A remark on computing distance
*         functions",J.Comp.Phys.,Vol.163,2000,pp.51-67
*           input    :sca, itnum
*           output   :dist
*           temporary:stmp
*           flagging :iflag
*******************************************************************************/
  boil::timer.start("cipcsl2 distfunc");

  //nlayer= 12; /* Nr. layer to be calculated around the free-surface */
  //nlayer= 16; /* Nr. layer to be calculated around the free-surface */
  const real phimin1=-dxmin * nlayer;
  const real phimax1= dxmin * nlayer;
  const int ifmax=nlayer+4;
  const int ifmin=-nlayer-4;

  /*-------------+
  |  initialize  |
  +-------------*/
  for_aijk(i,j,k){
    fn[i][j][k]=0.0;
  }

  /*-----------+
  |  set flag  |
  +-----------*/
  set_iflag();

#if 0
  boil::plot->plot(clr,iflag,"clr-iflag", time->current_step());
  exit(0);
#endif

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

  for(int it=1; it<=itnum; it++){
#if 0
    std::cout<<"dist: "<<it<<" "<<dist[4][17][1]<<" "<<sca[4][17][1]<<"\n";
    std::cout<<"dflag= "<<iflag[4][17][1]<<" "<<nlayer<<"\n";
#endif
    for_ijk(i,j,k){
      if(abs(iflag[i][j][k])<=nlayer){  // only near region to interface
        real dx = boil::minr(phi.dxc(i),phi.dyc(j),phi.dzc(k)); // crude
        real dtau = 0.5*dx;
        real asgn=copysign(1.0, stmp[i][j][k]);
        if(iflag[i][j][k]==0){
          /* cells next to free surface (NFCell) */
          /* x-direction */
          real dpdx=0.0;
          if(iflag[i-1][j][k]==0 && iflag[i+1][j][k]==0) {
            dpdx=(stmp[i+1][j][k]-stmp[i-1][j][k])/(phi.dxw(i)+phi.dxe(i));
          } else if(iflag[i-1][j][k]==0){
            dpdx=(stmp[i][j][k]-stmp[i-1][j][k])/phi.dxw(i);
          } else if (iflag[i+1][j][k]==0) {
            dpdx=(stmp[i+1][j][k]-stmp[i][j][k])/phi.dxe(i);
          }
          /* y-direction */
          real dpdy=0.0;
          if(iflag[i][j-1][k]==0 && iflag[i][j+1][k]==0) {
            dpdy=(stmp[i][j+1][k]-stmp[i][j-1][k])/(phi.dys(j)+phi.dyn(j));
          } else if(iflag[i][j-1][k]==0){
            dpdy=(stmp[i][j][k]-stmp[i][j-1][k])/phi.dys(j);
          } else if (iflag[i][j+1][k]==0) {
            dpdy=(stmp[i][j+1][k]-stmp[i][j][k])/phi.dyn(j);
          }
          /* z-direction */
          real dpdz=0.0;
          if(iflag[i][j][k-1]==0 && iflag[i][j][k+1]==0) {
            dpdz=(stmp[i][j][k+1]-stmp[i][j][k-1])/(phi.dzb(k)+phi.dzt(k));
          } else if(iflag[i][j][k-1]==0){
            dpdz=(stmp[i][j][k]-stmp[i][j][k-1])/phi.dzb(k);
          } else if (iflag[i][j][k+1]==0) {
            dpdz=(stmp[i][j][k+1]-stmp[i][j][k])/phi.dzt(k);
          }

          real dpdxyz=sqrt(dpdx*dpdx+dpdy*dpdy+dpdz*dpdz);
          dpdxyz = std::max(dpdxyz,1.0e-12);
          real d = stmp[i][j][k] / dpdxyz;
          fn[i][j][k] = dtau /dx * ( asgn * fabs(dist[i][j][k]) - d);
#if 0
          if(i==4&&j==17&&k==1) {
            std::cout<<"dtam= "<<dtau<<" "<<dx<<" "<<asgn<<" "
                     <<dist[i][j][k]<<" "<<stmp[i][j][k]<<" "<<dpdxyz<<"\n";
            std::cout<<"dpdx: "<<dpdx<<" "<<dpdy<<" "<<dpdz<<"\n";
            std::cout<<"stmp: "<<stmp[i][j][k]<<" "<<stmp[i][j][k-1]<<"\n";
            std::cout<<"iflag: "<<iflag[i][j][k-1]<<" "<<iflag[i][j][k+1]<<"\n";
            exit(0);
          }
#endif
        } else {
          /* cells except NFCell */
          real a=(dist[i][j][k]-dist[i-1][j][k])/phi.dxw(i);
          real b=(dist[i+1][j][k]-dist[i][j][k])/phi.dxe(i);
          real c=(dist[i][j][k]-dist[i][j-1][k])/phi.dys(j);
          real d=(dist[i][j+1][k]-dist[i][j][k])/phi.dyn(j);
          real e=(dist[i][j][k]-dist[i][j][k-1])/phi.dzb(k);
          real f=(dist[i][j][k+1]-dist[i][j][k])/phi.dzt(k);
          real g;

	  if(iflag[i-1][j][k]<-1000)a=0;
	  if(iflag[i+1][j][k]<-1000)b=0;
	  if(iflag[i][j-1][k]<-1000)c=0;
	  if(iflag[i][j+1][k]<-1000)d=0;
	  if(iflag[i][j][k-1]<-1000)e=0;
	  if(iflag[i][j][k+1]<-1000)f=0;

	  if(iflag[i][j][k]<-1000){
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
          fn[i][j][k] = dtau * asgn * g;
        }
      }
    }

#if 0
    std::cout<<"update= "<<fn[4][17][1]<<"\n";
#endif
    /* update Jakobi */
    for_ijk(i,j,k){
      if(abs(iflag[i][j][k])<=nlayer){
        dist[i][j][k] -= fn[i][j][k];
      }
    }

    insert_bc_dist(dist);
    dist.exchange_all();
  }

  for_aijk(i,j,k){
    if(int(iflag[i][j][k])==ifmax){
      dist[i][j][k]= real(ifmax)*dxmin;
    } else if(int(iflag[i][j][k])==-ifmax){
      dist[i][j][k]= real(ifmin)*dxmin;
    }
  }
  insert_bc_dist(dist);
  dist.exchange_all();

#if 0
  boil::plot->plot(clr,dist,iflag,
    "distfunc_clr-dist-iflag", time->current_step());
  exit(0);
#endif
  boil::timer.stop("cipcsl2 distfunc");

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_distfunc.cpp,v 1.3 2015/05/05 15:19:13 sato Exp $'/
+-----------------------------------------------------------------------------*/
