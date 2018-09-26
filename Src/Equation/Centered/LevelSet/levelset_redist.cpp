#include "levelset.h"
#include <iomanip>

/******************************************************************************/
void LevelSet::redist(const int itnum) {
/***************************************************************************//**
*  /brief Convert color function to distance function.
*         Reference:G.,Russo and P.,Smereka,"A remark on computing distance
*         functions",J.Comp.Phys.,Vol.163,2000,pp.51-67
*           input             :itnum
*           distance function :phi
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
    kappa[i][j][k]=0.0;
  }

  /*-----------+
  |  set flag  |
  +-----------*/
  for_aijk(i,j,k) {
    if(phi[i][j][k]<phisurf){
      dflag[i][j][k]=ifmin;
    } else {
      dflag[i][j][k]=ifmax;
    }
  }

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=0; i<ni()-1; i++){
    for_jk(j,k){
      if((phi[i][j][k]-phisurf)*(phi[i+1][j][k]-phisurf)<=0.0){
         dflag[i  ][j][k]=0;
         dflag[i+1][j][k]=0;
      }
    }
  }
  /* j-direction */
  for(int j=0; j<nj()-1; j++){
    for_ik(i,k){
      if((phi[i][j][k]-phisurf)*(phi[i][j+1][k]-phisurf)<=0.0){
        dflag[i][j  ][k]=0;
        dflag[i][j+1][k]=0;
      }
    }
  }
  /* k-direction */
  for(int k=0; k<nk()-1; k++){
    for_ij(i,j){
      if((phi[i][j][k]-phisurf)*(phi[i][j][k+1]-phisurf)<=0.0){
        dflag[i][j][k  ]=0;
        dflag[i][j][k+1]=0;
      }
    }
  }
  insert_bc_dist(dflag);
  dflag.exchange();

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
        if(fabs(dflag[i][j][k])==ifmax&&fabs(dflag[i+1][j][k])==(layer-1)){
           dflag[i  ][j][k]=layer*int(copysign(1.0,dflag[i  ][j][k]));
        } else if(fabs(dflag[i+1][j][k])==ifmax&&fabs(dflag[i][j][k])==(layer-1)){
           dflag[i+1][j][k]=layer*int(copysign(1.0,dflag[i+1][j][k]));
        }
      }
    }
    /* j-direction */
    for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
        if(fabs(dflag[i][j][k])==ifmax&&fabs(dflag[i][j+1][k])==(layer-1)){
           dflag[i  ][j][k]=layer*int(copysign(1.0,dflag[i  ][j][k]));
        } else if(fabs(dflag[i][j+1][k])==ifmax&&fabs(dflag[i][j][k])==(layer-1)){
           dflag[i][j+1][k]=layer*int(copysign(1.0,dflag[i][j+1][k]));
        }
      }
    }
    /* k-direction */
    for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
        if(fabs(dflag[i][j][k])==ifmax&&fabs(dflag[i][j][k+1])==(layer-1)){
           dflag[i  ][j][k]=layer*int(copysign(1.0,dflag[i  ][j][k]));
        } else if(fabs(dflag[i][j][k+1])==ifmax&&fabs(dflag[i][j][k])==(layer-1)){
           dflag[i][j][k+1]=layer*int(copysign(1.0,dflag[i][j][k+1]));
        }
      }
    }
    insert_bc_dist(dflag);
    dflag.exchange();
  }
  //boil::plot->plot(clr,dist,dflag, "clr-dist-dflag", time->current_step());

  /* store initial distance function */
  for_aijk(i,j,k){
    stmp[i][j][k]=phi[i][j][k];
  }

  for(int it=1; it<=itnum; it++){
    /* normal loop */
    for_ijk(i,j,k){
      if(fabs(dflag[i][j][k])<=nlayer){
        real dx = std::min(phi.dzc(k),std::min(phi.dxc(i),phi.dyc(j))); // crude
        real dtau = 0.5*dx;
        real asgn=copysign(1.0, stmp[i][j][k]);
        if(dflag[i][j][k]==0){
          /* cells next to free surface */
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
          kappa[i][j][k] = dtau /dx * ( asgn * fabs(phi[i][j][k]) - d);
        } else {
          /* cells except NFCell */
          real a=(phi[i][j][k]-phi[i-1][j][k])/phi.dxw(i);
          real b=(phi[i+1][j][k]-phi[i][j][k])/phi.dxe(i);
          real c=(phi[i][j][k]-phi[i][j-1][k])/phi.dys(j);
          real d=(phi[i][j+1][k]-phi[i][j][k])/phi.dyn(j);
          real e=(phi[i][j][k]-phi[i][j][k-1])/phi.dzb(k);
          real f=(phi[i][j][k+1]-phi[i][j][k])/phi.dzt(k);
          real g;
          if(phi[i][j][k]>=0){
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
          kappa[i][j][k] = dtau * asgn * g;
        }
      }
    }

    /* update Gauss-Saidel */
    for_ijk(i,j,k){
      if(fabs(dflag[i][j][k])<=nlayer){
        phi[i][j][k] -= kappa[i][j][k];
      }
    }

    insert_bc_dist(phi);
    phi.exchange_all();
  }

  for_aijk(i,j,k){
    if(dflag[i][j][k]==ifmax){
      phi[i][j][k]= real(ifmax)*dxmin;
    } else if(dflag[i][j][k]==-ifmax){
      phi[i][j][k]= real(ifmin)*dxmin;
    }
  }
  insert_bc_dist(phi);
  phi.exchange_all();

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: levelset_redist.cpp,v 1.2 2012/09/13 08:42:27 niceno Exp $'/
+-----------------------------------------------------------------------------*/
