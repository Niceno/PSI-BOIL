#include "pressure.h"

real kappa_avei(const Scalar & ka, int i, int j, int k);
real kappa_avej(const Scalar & ka, int i, int j, int k);
real kappa_avek(const Scalar & ka, int i, int j, int k);

/***************************************************************************//**
*  Calculate right hand side for ghost fluid method
*******************************************************************************/
void Pressure::ghost(const Scalar & c, const Scalar & kappa) {

  real sigma = fluid()->sigma()->value();
  real cint = 0.5;

  if( dom->ibody().nccells() == 0 ) {
    for_ijk(i,j,k) {

      real cc = c[i][j][k];
      real cw = c[i-1][j][k];
      real ce = c[i+1][j][k];
      real cs = c[i][j-1][k];
      real cn = c[i][j+1][k];
      real cb = c[i][j][k-1];
      real ct = c[i][j][k+1];

      /* west */
      if ((cc-cint)*(cw-cint)<0.0) {
        //std::cout<<"ghost:w";
        Comp m =  Comp::u();
        real rho = fluid()->rho(m,i,j,k);
        real s   = dSx(Sign::neg(),i,j,k);
        real dd  = dxw(i);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avei(kappa,i,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* east */
      if ((cc-cint)*(ce-cint)<0.0) {
        Comp m =  Comp::u();
        real rho = fluid()->rho(m,i+1,j,k);
        real s   = dSx(Sign::neg(),i+1,j,k);
        real dd  = dxe(i);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avei(kappa,i+1,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* south */
      if ((cc-cint)*(cs-cint)<0.0) {
        Comp m =  Comp::v();
        real rho = fluid()->rho(m,i,j,k);
        real s   = dSy(Sign::neg(),i,j,k);
        real dd  = dys(j);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avej(kappa,i,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* north */
      if ((cc-cint)*(cn-cint)<0.0) {
        Comp m =  Comp::v();
        real rho = fluid()->rho(m,i,j+1,k);
        real s   = dSy(Sign::neg(),i,j+1,k);
        real dd  = dys(j);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avej(kappa,i,j+1,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* bottom */
      if ((cc-cint)*(cb-cint)<0.0) {
        Comp m =  Comp::w();
        real rho = fluid()->rho(m,i,j,k);
        real s   = dSz(Sign::neg(),i,j,k);
        real dd  = dzb(k);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avek(kappa,i,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* top */
      if ((cc-cint)*(ct-cint)<0.0) {
        Comp m =  Comp::w();
        real rho = fluid()->rho(m,i,j,k+1);
        real s   = dSz(Sign::neg(),i,j,k+1);
        real dd  = dzt(k);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avek(kappa,i,j,k+1);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }
    }
  } else {

    for_ijk(i,j,k) {

      if( dom->ibody().off_p(i,j,k)) continue;

      real cc = c[i][j][k];
      real cw = c[i-1][j][k];
      real ce = c[i+1][j][k];
      real cs = c[i][j-1][k];
      real cn = c[i][j+1][k];
      real cb = c[i][j][k-1];
      real ct = c[i][j][k+1];

      /* west */
      if ((cc-cint)*(cw-cint)<0.0) {
        //std::cout<<"ghost:w";
        Comp m =  Comp::u();
        real rho = fluid()->rho(m,i,j,k);
        real s   = dSx(Sign::neg(),i,j,k);
        real dd  = dxw(i);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avei(kappa,i,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

        /* east */
      if ((cc-cint)*(ce-cint)<0.0) {
        Comp m =  Comp::u();
        real rho = fluid()->rho(m,i+1,j,k);
        real s   = dSx(Sign::neg(),i+1,j,k);
        real dd  = dxe(i);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avei(kappa,i+1,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* south */
      if ((cc-cint)*(cs-cint)<0.0) {
        Comp m =  Comp::v();
        real rho = fluid()->rho(m,i,j,k);
        real s   = dSy(Sign::neg(),i,j,k);
        real dd  = dys(j);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avej(kappa,i,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* north */
      if ((cc-cint)*(cn-cint)<0.0) {
        Comp m =  Comp::v();
        real rho = fluid()->rho(m,i,j+1,k);
        real s   = dSy(Sign::neg(),i,j+1,k);
        real dd  = dys(j);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avej(kappa,i,j+1,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }

      /* bottom */
      if ((cc-cint)*(cb-cint)<0.0) {
        Comp m =  Comp::w();
        real rho = fluid()->rho(m,i,j,k);
        real s   = dSz(Sign::neg(),i,j,k);
        real dd  = dzb(k);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avek(kappa,i,j,k);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }
  
      /* top */
      if ((cc-cint)*(ct-cint)<0.0) {
        Comp m =  Comp::w();
        real rho = fluid()->rho(m,i,j,k+1);
        real s   = dSz(Sign::neg(),i,j,k+1);
        real dd  = dzt(k);
        int iphase = 1;
        if (cc<cint) iphase = -1;
        real kap = kappa_avek(kappa,i,j,k+1);
        fext[i][j][k] += real(iphase)*s*sigma*kap/(rho*dd);
      }
    }
  }
}

/* harmonic averaging --------------------------------------------------------*/
real kappa_avei(const Scalar & kappa, int i, int j, int k){
  real kappa_ave = 0.5*(kappa[i-1][j][k]+kappa[i][j][k]);
  if (kappa[i-1][j][k]*kappa[i][j][k]>0.0) {
    kappa_ave = 2.0 * kappa[i-1][j][k] * kappa[i][j][k]
                   / (kappa[i-1][j][k] + kappa[i][j][k]);
  }
  return kappa_ave;
}
real kappa_avej(const Scalar & kappa, int i, int j, int k){
  real kappa_ave = 0.5*(kappa[i][j-1][k]+kappa[i][j][k]);
  if (kappa[i][j-1][k]*kappa[i][j][k]>0.0) {
    kappa_ave = 2.0 * kappa[i][j-1][k] * kappa[i][j][k]
                   / (kappa[i][j-1][k] + kappa[i][j][k]);
  }
  return kappa_ave;
}
real kappa_avek(const Scalar & kappa, int i, int j, int k){
  real kappa_ave = 0.5*(kappa[i][j][k-1]+kappa[i][j][k]);
  if (kappa[i][j][k-1]*kappa[i][j][k]>0.0) {
    kappa_ave = 2.0 * kappa[i][j][k-1] * kappa[i][j][k]
                   / (kappa[i][j][k-1] + kappa[i][j][k]);
  }
  return kappa_ave;
}
