#include "phasechange.h"
using namespace std;

/******************************************************************************/
void PhaseChange::gradtx5( const int i, const int j, const int k,
              real * txm, real * txp, const int mm) const {
/***************************************************************************//**
*  \brief calculate grad(tpr) in x-direction
*         mm=-1: grad(tpr) between m and interface
*         mm= 1: grad(tpr) between interface and p
*           m: [i  ][j][k]
*           p: [i+1][j][k]
*******************************************************************************/
  real clrm=std::max(0.0,std::min(1.0,clr[i][j][k]));
  real clrp=std::max(0.0,std::min(1.0,clr[i+1][j][k]));

  real wm = (phisurf-clrm)/(clrp-clrm);
  real wp = 1.0-wm;

  real dxm, dxp, tprm, tprp;
  if(mm<0){
    if (wm>=epsl) {
      real a = clr.dxw(i);
      real b = clr.dxe(i) * wm;
#if 0
      * txm = (-b*b*tpr[i-1][j][k]
               +(b*b-a*a)*tpr[i][j][k]
               +a*a*tsat)
              /(a*b*(a+b));
#else
      * txm = b*b*(tpr[i][j][k]-tpr[i-1][j][k])
             +a*a*(tsat-tpr[i][j][k]);
      * txm /= (a*b*(a+b));
#endif
    } else {
      dxm = clr.dxw(i) + clr.dxe(i) * wm;
      tprm = tpr[i-1][j][k];
      * txm = (tsat-tprm)/dxm;
    }
  } else {
    if (wp>=epsl) {
      real a = clr.dxe(i) * wp;
      real b = clr.dxe(i+1);
#if 0
      * txp = (-b*b*tsat
               +(b*b-a*a)*tpr[i+1][j][k]
               +a*a*tpr[i+2][j][k])
              /(a*b*(a+b));
#else
      * txp = b*b*(tpr[i+1][j][k]-tsat)
             +a*a*(tpr[i+2][j][k]-tpr[i+1][j][k]);
      * txp /= (a*b*(a+b));
#endif
    } else {
      dxp = clr.dxe(i+1) + clr.dxe(i) * wp;
      tprp = tpr[i+2][j][k];
      * txp = (tprp-tsat)/dxp;
    }
  }
}

/******************************************************************************/
void PhaseChange::gradty5( const int i, const int j, const int k,
              real * tym, real * typ, const int mm) const {
/***************************************************************************//**
*  \brief calculate grad(tpr) in y-direction
*         mm=-1: grad(tpr) between m and interface
*         mm= 1: grad(tpr) between interface and p
*           m: [i][j  ][k]
*           p: [i][j+1][k]
*******************************************************************************/
  real clrm=std::max(0.0,std::min(1.0,clr[i][j][k]));
  real clrp=std::max(0.0,std::min(1.0,clr[i][j+1][k]));

  real wm = (phisurf-clrm)/(clrp-clrm);
  real wp = 1.0-wm;

  real dym, dyp, tprm, tprp;
  if(mm<0){
    if(wm>epsl){
      real a = clr.dys(j);
      real b = clr.dyn(j) * wm;
#if 0
      * tym = (-b*b*tpr[i][j-1][k]
               +(b*b-a*a)*tpr[i][j][k]
               +a*a*tsat)
              /(a*b*(a+b));
#else
      * tym = b*b*(tpr[i][j][k]-tpr[i][j-1][k])
             +a*a*(tsat-tpr[i][j][k]);
      * tym /= (a*b*(a+b));
#endif
    } else {
      dym = clr.dys(j) + clr.dyn(j) * wm;
      tprm = tpr[i][j-1][k];
      * tym = (tsat-tprm)/dym;
    }
  } else {
    if(wp>epsl){
      real a = clr.dyn(j) * wp;
      real b = clr.dyn(j+1);
#if 0
      * typ = (-b*b*tsat
               +(b*b-a*a)*tpr[i][j+1][k]
               +a*a*tpr[i][j+2][k])
              /(a*b*(a+b));
#else
      * typ = b*b*(tpr[i][j+1][k]-tsat)
             +a*a*(tpr[i][j+2][k]-tpr[i][j+1][k]);
      * typ /= (a*b*(a+b));
#endif
    } else {
      dyp = clr.dyn(j+1) + clr.dyn(j) * wp;
      tprp = tpr[i][j+2][k];
      * typ = (tprp-tsat)/dyp;
    }
  }
}

/******************************************************************************/
void PhaseChange::gradtz5( const int i, const int j, const int k,
              real * tzm, real * tzp, const int mm) const {
/***************************************************************************//**
*  \brief calculate grad(tpr) in z-direction
*         mm=-1: grad(tpr) between m and interface
*         mm= 1: grad(tpr) between interface and p
*           m: [i][j][k  ]
*           p: [i][j][k+1]
*******************************************************************************/
  real clrm=std::max(0.0,std::min(1.0,clr[i][j][k]));
  real clrp=std::max(0.0,std::min(1.0,clr[i][j][k+1]));

  real wm = (phisurf-clrm)/(clrp-clrm);
  real wp = 1.0-wm;

  real dzm, dzp, tprm, tprp;
  if(mm<0){
    if(wm>epsl){
      real a = clr.dzb(k);
      real b = clr.dzt(k) * wm;
#if 0
      * tzm = (-b*b*tpr[i][j][k-1]
               +(b*b-a*a)*tpr[i][j][k]
               +a*a*tsat)
              /(a*b*(a+b));
#else
      * tzm = b*b*(tpr[i][j][k]-tpr[i][j][k-1])
             +a*a*(tsat-tpr[i][j][k]);
      * tzm /= (a*b*(a+b));
#endif
    } else {
      dzm = clr.dzb(k) + clr.dzt(k) * wm;
      tprm = tpr[i][j][k-1];
      * tzm = (tsat-tprm)/dzm;
    }
  } else {
    if(wp>epsl){
      real a = clr.dzt(k) * wp;
      real b = clr.dzt(k+1);
#if 0
      * tzp = (-b*b*tsat
               +(b*b-a*a)*tpr[i][j][k+1]
               +a*a*tpr[i][j][k+2])
              /(a*b*(a+b));
#else
      * tzp = b*b*(tpr[i][j][k+1]-tsat)
             +a*a*(tpr[i][j][k+2]-tpr[i][j][k+1]);
      * tzp /= (a*b*(a+b));
#endif
    } else {
      dzp = clr.dzt(k+1) + clr.dzt(k) * wp;
      tprp = tpr[i][j][k+2];
      * tzp = (tprp-tsat)/dzp;
    }
  } 
}
