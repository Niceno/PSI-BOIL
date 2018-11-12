#include "cipcsl2.h"
#include <iomanip>
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::tension(Vector * vec, const Matter matt) {
/******************************************************************************/
  tension(vec, matt, phi);
}
/******************************************************************************/
void CIPCSL2::tension(Vector * vec, const Matter matt, Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate surface tension
*         Algorithm
*           1st step: calculate distance or smoothed color function
*           2nd step: calculate curvature from the distance/smoothed function
*           3rd step: calculate body force
*         Variables
*           color function          : clr
*           smoothed color function : sclr 
*           distance function       : dist
*           curvature               : kappa
*           body force              : vec
*******************************************************************************/
  boil::timer.start("cipcsl2 tension");
#ifdef DEBUG
  std::cout<<"tension::begin "<<boil::cart.iam()<<"\n";
#endif

  /*-------------------+
  |  1st and 2nd step  |
  +-------------------*/
  curvature();

  /*-----------+
  |  3rd step  |
  +-----------*/
  real rho_diff = matt.rho(1)-matt.rho(0);
  real rho_ave = 0.5*(matt.rho(1)+matt.rho(0));

  /* density-scaled distance function */
  if(use_dist_for_kappa){
    //real eps=1.5*dxmin;
    for_aijk(i,j,k){
      real eps=1.5*max(max(phi.dxc(i),phi.dyc(j)),phi.dzc(k));
      if(dist[i][j][k]<-eps){
        stmp[i][j][k]=0.0;
      } else if(dist[i][j][k]>eps) {
        stmp[i][j][k]=1.0;
      } else {
        stmp[i][j][k]= 0.5 + dist[i][j][k]/(2.0*eps)
                     + 1.0/(2.0*pi)*sin(pi*dist[i][j][k]/eps);
      }
    }
    for_aijk(i,j,k){
      stmp[i][j][k] = matt.rho(0)*(1.0-stmp[i][j][k])
                    + matt.rho(1)*stmp[i][j][k];
    }
  } else {
    for_aijk(i,j,k){
      stmp[i][j][k] = matt.rho(i,j,k);
    }
  }

  /* calculate body force */
  Comp m;
  m = Comp::u();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        real kappa_ave = 0.5*(kappa[i-1][j][k]+kappa[i][j][k]);
        if (kappa[i-1][j][k]*kappa[i][j][k]>0.0) {
          kappa_ave = 2.0 * kappa[i-1][j][k] * kappa[i][j][k]
                       / (kappa[i-1][j][k] + kappa[i][j][k]);
        }
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave
              * (sca[i][j][k] - sca[i-1][j][k])/vec->dxc(m,i)
              * vec->dV(m,i,j,k);
      }
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        real kappa_ave = 0.5*(kappa[i-1][j][k]+kappa[i][j][k]);
        if (kappa[i-1][j][k]*kappa[i][j][k]>0.0) {
          kappa_ave = 2.0 * kappa[i-1][j][k] * kappa[i][j][k]
                       / (kappa[i-1][j][k] + kappa[i][j][k]);
        }
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave
              * (stmp[i][j][k]-stmp[i-1][j][k])/vec->dxc(m,i)
              / rho_diff * 0.5*(stmp[i][j][k]+stmp[i-1][j][k])
              / rho_ave
              * vec->dV(m,i,j,k);
      }
    }
  }

  m = Comp::v();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        real kappa_ave = 0.5*(kappa[i][j-1][k]+kappa[i][j][k]);
        if (kappa[i][j-1][k]*kappa[i][j][k]>0.0) {
          kappa_ave = 2.0 * kappa[i][j-1][k] * kappa[i][j][k]
                       / (kappa[i][j-1][k] + kappa[i][j][k]);
        }
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave
              * (sca[i][j][k] - sca[i][j-1][k])/vec->dyc(m,j)
              * vec->dV(m,i,j,k);
      }
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        real kappa_ave = 0.5*(kappa[i][j-1][k]+kappa[i][j][k]);
        if (kappa[i][j-1][k]*kappa[i][j][k]>0.0) {
          kappa_ave = 2.0 * kappa[i][j-1][k] * kappa[i][j][k]
                       / (kappa[i][j-1][k] + kappa[i][j][k]);
        }
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave
              * (stmp[i][j][k]-stmp[i][j-1][k])/vec->dyc(m,j)
              / rho_diff * 0.5*(stmp[i][j][k]+stmp[i][j-1][k])
              / rho_ave
              * vec->dV(m,i,j,k);
      }
    }
  }

  m = Comp::w();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        real kappa_ave = 0.5*(kappa[i][j][k-1]+kappa[i][j][k]);
        if (kappa[i][j][k-1]*kappa[i][j][k]>0.0) {
          kappa_ave = 2.0 * kappa[i][j][k-1] * kappa[i][j][k]
                       / (kappa[i][j][k-1] + kappa[i][j][k]);
        }
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave
              * (sca[i][j][k] - sca[i][j][k-1])/vec->dzc(m,k)
              * vec->dV(m,i,j,k);
      }
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        real kappa_ave = 0.5*(kappa[i][j][k-1]+kappa[i][j][k]);
        if (kappa[i][j][k-1]*kappa[i][j][k]>0.0) {
          kappa_ave = 2.0 * kappa[i][j][k-1] * kappa[i][j][k]
                       / (kappa[i][j][k-1] + kappa[i][j][k]);
        }
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave
              * (stmp[i][j][k]-stmp[i][j][k-1])/vec->dzc(m,k)
              / rho_diff * 0.5*(stmp[i][j][k]+stmp[i][j][k-1])
              / rho_ave
              * vec->dV(m,i,j,k);
      }
    }
  }
  vec->exchange();

  boil::timer.stop("cipcsl2 tension");
}
