#include "cipcsl2.h"
#include <iomanip>
//#define ORIGINAL_CSF
#define KAPPA_DIFFEQ

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
*           smoothed color function : stmp
*           distance function       : dist
*           curvature               : kappa
*           body force              : vec
*******************************************************************************/
  boil::timer.start("cipcsl2 tension");

  /*-----------+
  |  1st step  |
  +-----------*/
#ifdef KAPPA_DIFFEQ
  smooth(clr, stmp, itsmooth);
  //for_aijk(i,j,k)stmp[i][j][k]=clr[i][j][k];
#else
  distfunc(clr, 24);
#endif

  /*-----------+
  |  2nd step  |
  +-----------*/
  /* calculate curvature */
#ifdef KAPPA_DIFFEQ
  curv(stmp);
#else
  curv(dist);
#endif

#if 0
  if(time->current_step()==1||time->current_step()==100){
    boil::plot->plot(kappa,nx,ny,nz, "kappa-nx-ny-nz", time->current_step());
    boil::plot->plot(clr,kappa,stmp, "clr-kappa-stmp", time->current_step());
    boil::plot->plot(clr,kappa,dist, "clr-kappa-dist", time->current_step());
    exit(0);
  }
#endif

  /* wall adhesion */
#ifdef KAPPA_DIFFEQ
  bdcurv(stmp);
#else
  bdcurv(dist);
#endif
  insert_bc_kappa2(kappa);
  kappa.exchange();

#if 0
  if(time->current_step()==3||time->current_step()==100){
    boil::plot->plot(kappa,nx,ny,nz, "kappa-nx-ny-nz", time->current_step());
    boil::plot->plot(clr,kappa,stmp, "clr-kappa-stmp", time->current_step());
    boil::plot->plot(clr,dflag,stmp, "clr-dflag-stmp", time->current_step());
    exit(0);
  }
#endif
  //boil::plot->plot(clr,kappa,dist, 
  //  "tension_clr-kappa-dist", time->current_step());
  //exit(0);

  /*-----------+
  |  3rd step  |
  +-----------*/
  real rho_diff = fabs(matt.rho(1)-matt.rho(0));
  real rho_ave = 0.5*(matt.rho(1)+matt.rho(0));

#ifndef ORIGINAL_CSF
  if(i_st_dist==1){
    // use density based on distance function
    real eps=eps_st*dxmin;
    for_aijk(i,j,k){
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
    //boil::plot->plot(clr,dist,stmp, "clr-dist-stmp", time->current_step());
  } else {
    // use density based on color function
    for_aijk(i,j,k){
      stmp[i][j][k] = matt.rho(i,j,k);
    }
  }
#endif

  /* calculate body force */
  Comp m;
  m = Comp::u();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
#ifdef IB
      if(dom->ibody().on(m,i,j,k))
#endif
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i-1][j][k]+kappa[i][j][k])
                         * (sca[i][j][k] - sca[i-1][j][k])/vec->dxc(m,i)
                         * vec->dV(m,i,j,k);
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
#ifdef IB
      if(dom->ibody().on(m,i,j,k))
#endif
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i-1][j][k]+kappa[i][j][k])
#ifdef ORIGINAL_CSF
                         * (sca[i][j][k] - sca[i-1][j][k])/vec->dxc(m,i)
#else
                         * (stmp[i][j][k]-stmp[i-1][j][k])/vec->dxc(m,i)
                         /// rho_diff * matt.rho(m,i,j,k) / rho_ave
                         / rho_diff * 0.5*(stmp[i][j][k]+stmp[i-1][j][k])
                         / rho_ave
#endif
                       * vec->dV(m,i,j,k);
    }
  }

  m = Comp::v();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
#ifdef IB
      if(dom->ibody().on(m,i,j,k))
#endif
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j-1][k]+kappa[i][j][k])
                         * (sca[i][j][k] - sca[i][j-1][k])/vec->dyc(m,j)
                         * vec->dV(m,i,j,k);
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
#ifdef IB
      if(dom->ibody().on(m,i,j,k))
#endif
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j-1][k]+kappa[i][j][k])
#ifdef ORIGINAL_CSF
                         * (sca[i][j][k] - sca[i][j-1][k])/vec->dyc(m,j)
#else
                         * (stmp[i][j][k]-stmp[i][j-1][k])/vec->dyc(m,j)
                         /// rho_diff * matt.rho(m,i,j,k) / rho_ave
                         / rho_diff * 0.5*(stmp[i][j][k]+stmp[i][j-1][k])
                         / rho_ave
#endif
                         * vec->dV(m,i,j,k);
    }
  }

  m = Comp::w();
  if(rho_diff==0.0){
    for_vmijk((*vec),m,i,j,k) {
#ifdef IB
      if(dom->ibody().on(m,i,j,k))
#endif
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j][k-1]+kappa[i][j][k])
                         * (sca[i][j][k] - sca[i][j][k-1])/vec->dzc(m,k)
                         * vec->dV(m,i,j,k);
    }
  } else {
    for_vmijk((*vec),m,i,j,k) {
#ifdef IB
      if(dom->ibody().on(m,i,j,k))
#endif
      (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                         * 0.5*(kappa[i][j][k-1]+kappa[i][j][k])
#ifdef ORIGINAL_CSF
                         * (sca[i][j][k] - sca[i][j][k-1])/vec->dzc(m,k)
#else 
                         * (stmp[i][j][k]-stmp[i][j][k-1])/vec->dzc(m,k)
                         /// rho_diff * matt.rho(m,i,j,k) / rho_ave
                         / rho_diff * 0.5*(stmp[i][j][k]+stmp[i][j][k-1])
                         / rho_ave
#endif
                         * vec->dV(m,i,j,k);
    }
  }
  vec->exchange();

  boil::timer.stop("cipcsl2 tension");
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_tension.cpp,v 1.1 2015/02/19 09:41:42 sato Exp $'/
+-----------------------------------------------------------------------------*/
