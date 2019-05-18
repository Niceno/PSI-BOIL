#include "vof.h"
#include <iomanip>

/******************************************************************************/
void VOF::tension(Vector * vec, const Matter matt) {
/******************************************************************************/

  /*----------------------------------+
  |  1st step: curvature calculation  |
  +----------------------------------*/
  curvature();
#if 0
  diffuse(kap,stmp4,4);
  for_avijk(kap,i,j,k) {
    //stmp2[i][j][k] = kap[i][j][k];
    kap[i][j][k] = stmp4[i][j][k];
  }

  for_avijk(kap,i,j,k) {
    //kap[i][j][k] = 1./1e-3;
  }
#endif

  tension(vec,matt,kappa);
}

/******************************************************************************/
void VOF::tension(Vector * vec, const Matter matt, const Scalar & kap) {
/***************************************************************************//**
*  \brief Calculate surface tension
*         Algorithm
*           1st step: calculate curvature
*           2nd step: calculate body force
*         Variables
*           color function          : phi
*           curvature               : kap
*           body force              : vec
*******************************************************************************/
  boil::timer.start("vof tension");

  Comp m;

#if 0
  m = Comp::u();
  for_vmijk((*vec),m,i,j,k) {
    if(adens[i-1][j][k]>0.0||adens[i][j][k]>0.0) {
      real bdphi;
      if(bndclr) {
        bdphi = (*bndclr)[m][i][j][k];
      } else {
        bdphi = 0.5*(phi[i-1][j][k]+phi[i][j][k]);
      }
      real mult = bdphi*matt.sigma(m,i,j,k)*vec->dV(m,i,j,k);
      (*vec)[m][i][j][k] -= mult*(kap[i][j][k]-kap[i-1][j][k])/vec->dxc(m,i);
    }
  }

  m = Comp::v();
  for_vmijk((*vec),m,i,j,k) {
    if(adens[i][j-1][k]>0.0||adens[i][j][k]>0.0) {
      real bdphi;
      if(bndclr) {
        bdphi = (*bndclr)[m][i][j][k];
      } else {
        bdphi = 0.5*(phi[i][j-1][k]+phi[i][j][k]);
      }
      real mult = bdphi*matt.sigma(m,i,j,k)*vec->dV(m,i,j,k);
      (*vec)[m][i][j][k] -= mult*(kap[i][j][k]-kap[i][j-1][k])/vec->dyc(m,j);
    }
  }

  m = Comp::w();
  for_vmijk((*vec),m,i,j,k) {
    if(adens[i][j][k-1]>0.0||adens[i][j][k]>0.0) {
      real bdphi;
      if(bndclr) {
        bdphi = (*bndclr)[m][i][j][k];
      } else {
        bdphi = 0.5*(phi[i][j][k-1]+phi[i][j][k]);
      }
      real mult = bdphi*matt.sigma(m,i,j,k)*vec->dV(m,i,j,k);
      (*vec)[m][i][j][k] -= mult*(kap[i][j][k]-kap[i][j][k-1])/vec->dzc(m,k);
    }
  }

#else
  /*-----------------------+
  |  2nd step: body force  |
  +-----------------------*/
  real rho_diff = matt.rho(1)-matt.rho(0);
  real rho_ave = 0.5*(matt.rho(1)+matt.rho(0));

  if(rho_diff==0.0){
    m = Comp::u();
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave(kap[i-1][j][k],kap[i][j][k])
              * (phi[i][j][k] - phi[i-1][j][k])/vec->dxc(m,i)
              * vec->dV(m,i,j,k);
      }
    }
    m = Comp::v();
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave(kap[i][j-1][k],kap[i][j][k])
              * (phi[i][j][k] - phi[i][j-1][k])/vec->dyc(m,j)
              * vec->dV(m,i,j,k);
      }
    }
    m = Comp::w();
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave(kap[i][j][k-1],kap[i][j][k])
              * (phi[i][j][k] - phi[i][j][k-1])/vec->dzc(m,k)
              * vec->dV(m,i,j,k);
      }
    }
  } else {
    m = Comp::u();
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave(kap[i-1][j][k],kap[i][j][k])
              * (matt.rho(i,j,k)-matt.rho(i-1,j,k))/vec->dxc(m,i)
              / rho_diff * 0.5*(matt.rho(i,j,k)+matt.rho(i-1,j,k))
              / rho_ave
              * vec->dV(m,i,j,k);
      }
    }
    m = Comp::v();
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave(kap[i][j-1][k],kap[i][j][k])
              * (matt.rho(i,j,k)-matt.rho(i,j-1,k))/vec->dyc(m,j)
              / rho_diff * 0.5*(matt.rho(i,j,k)+matt.rho(i,j-1,k))
              / rho_ave
              * vec->dV(m,i,j,k);
      }
    }
    m = Comp::w();
    for_vmijk((*vec),m,i,j,k) {
      if(dom->ibody().on(m,i,j,k)) {
        (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
              * kappa_ave(kap[i][j][k-1],kap[i][j][k])
              * (matt.rho(i,j,k)-matt.rho(i,j,k-1))/vec->dzc(m,k)
              / rho_diff * 0.5*(matt.rho(i,j,k)+matt.rho(i,j,k-1))
              / rho_ave
              * vec->dV(m,i,j,k);
      }
    }
  }
#endif
  vec->exchange();

  boil::timer.stop("vof tension");
}
/******************************************************************************/
real VOF::kappa_ave(const real r1, const real r2) {
  real x;
  if       (!boil::realistic(r1)&&!boil::realistic(r2)) {
    x = 0.0;
  } else if(!boil::realistic(r1)) {
    x = r2;
  } else if(!boil::realistic(r2)) {
    x = r1;
  } else if(r1*r2>0.0) {
    x = 2.0 * r1 * r2 / (r1 + r2);
  } else {
    x = 0.5*(r1+r2);
  }
  return x;
}
