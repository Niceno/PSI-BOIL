#include "vof.h"
#include <iomanip>

void VOF::tension(Vector * vec, const Matter matt) {
  tension(vec,matt,phi);
}

/******************************************************************************/
void VOF::tension(Vector * vec, const Matter matt, const Scalar & scp) {
/***************************************************************************//**
*  \brief Calculate surface tension
*         Algorithm
*           1st step: calculate curvature
*           2nd step: calculate body force
*         Variables
*           color function          : scp
*           curvature               : kappa
*           body force              : vec
*******************************************************************************/
  boil::timer.start("vof tension");

  /*----------------------------------+
  |  1st step: curvature calculation  |
  +----------------------------------*/
  curvature();

  /*-----------------------+
  |  2nd step: body force  |
  +-----------------------*/
  real rho_diff = matt.rho(1)-matt.rho(0);
  real rho_ave = 0.5*(matt.rho(1)+matt.rho(0));

  Comp m;
  if(rho_diff==0.0){
    if(bflag_struct.ifull) {
      m = Comp::u();
      for_vmijk((*vec),m,i,j,k) {
        if(dom->ibody().on(m,i,j,k)) {
          (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                * kappa_ave(kappa[i-1][j][k],kappa[i][j][k])
                //* kappa_ave(kappa[i-1][j][k],kappa[i][j][k],tempflag[i-1][j][k],tempflag[i][j][k])
                * (scp[i][j][k] - scp[i-1][j][k])/vec->dxc(m,i)
                * vec->dV(m,i,j,k);
        }
      }
    }
    if(bflag_struct.jfull) {
      m = Comp::v();
      for_vmijk((*vec),m,i,j,k) {
        if(dom->ibody().on(m,i,j,k)) {
          (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                * kappa_ave(kappa[i][j-1][k],kappa[i][j][k])
                //* kappa_ave(kappa[i][j-1][k],kappa[i][j][k],tempflag[i][j-1][k],tempflag[i][j][k])
                * (scp[i][j][k] - scp[i][j-1][k])/vec->dyc(m,j)
                * vec->dV(m,i,j,k);
        }
      }
    }
    if(bflag_struct.kfull) {
      m = Comp::w();
      for_vmijk((*vec),m,i,j,k) {
        if(dom->ibody().on(m,i,j,k)) {
          (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                * kappa_ave(kappa[i][j][k-1],kappa[i][j][k])
                //* kappa_ave(kappa[i][j][k-1],kappa[i][j][k],tempflag[i][j][k-1],tempflag[i][j][k])
                * (scp[i][j][k] - scp[i][j][k-1])/vec->dzc(m,k)
                * vec->dV(m,i,j,k);
        }
      }
    }
  } else {
    if(bflag_struct.ifull) {
      m = Comp::u();
      for_vmijk((*vec),m,i,j,k) {
        if(dom->ibody().on(m,i,j,k)) {
          (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                * kappa_ave(kappa[i-1][j][k],kappa[i][j][k])
                //* kappa_ave(kappa[i-1][j][k],kappa[i][j][k],tempflag[i-1][j][k],tempflag[i][j][k])
                * (matt.rho(i,j,k)-matt.rho(i-1,j,k))/vec->dxc(m,i)
                / rho_diff * 0.5*(matt.rho(i,j,k)+matt.rho(i-1,j,k))
                / rho_ave
                * vec->dV(m,i,j,k);
        }
      }
    }
    if(bflag_struct.jfull) {
      m = Comp::v();
      for_vmijk((*vec),m,i,j,k) {
        if(dom->ibody().on(m,i,j,k)) {
          (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                * kappa_ave(kappa[i][j-1][k],kappa[i][j][k])
                //* kappa_ave(kappa[i][j-1][k],kappa[i][j][k],tempflag[i][j-1][k],tempflag[i][j][k])
                * (matt.rho(i,j,k)-matt.rho(i,j-1,k))/vec->dyc(m,j)
                / rho_diff * 0.5*(matt.rho(i,j,k)+matt.rho(i,j-1,k))
                / rho_ave
                * vec->dV(m,i,j,k);
        }
      }
    }
    if(bflag_struct.kfull) {
      m = Comp::w();
      for_vmijk((*vec),m,i,j,k) {
        if(dom->ibody().on(m,i,j,k)) {
          (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                * kappa_ave(kappa[i][j][k-1],kappa[i][j][k])
                //* kappa_ave(kappa[i][j][k-1],kappa[i][j][k],tempflag[i][j][k-1],tempflag[i][j][k])
                * (matt.rho(i,j,k)-matt.rho(i,j,k-1))/vec->dzc(m,k)
                / rho_diff * 0.5*(matt.rho(i,j,k)+matt.rho(i,j,k-1))
                / rho_ave
                * vec->dV(m,i,j,k);
        }
      }
    }
  }
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
/******************************************************************************/
real VOF::kappa_ave(const real r1, const real r2, const int i1, const int i2) {
  real x;
  //  be careful (Yohei) !!!
  if (i1==1 && i2==1) {
    if (r1*r2>0.0) {
      x = 2.0 * r1 * r2 / (r1 + r2);
    } else {
      x = 0.5*(r1+r2);
    }
  }
  if (i1==1 && i2==2) {
    x = r1;
  }
  if (i1==2 && i2==1) {
    x = r2;
  }
  if (i1==2 && i2==0) {
    x = r1;
  }
  if (i1==0 && i2==2) {
    x = r2;
  }
  if (i1==1 && i2==0) {
    std::cout<<"VOF::tension: error-kappa_ave\n";
    exit(0);
  }
  if (i1==0 && i2==1) {
    std::cout<<"VOF::tension: error-kappa_ave\n";
    exit(0);
  }

  return x;
}

