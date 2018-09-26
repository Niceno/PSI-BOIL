#include "cipcsl2.h"
#include "../../../Parallel/Out/out.h"
#include <iomanip>
using namespace std;

/******************************************************************************/
void CIPCSL2::curv_interface_ext() {
/***************************************************************************//**
*  \brief extrapolate kappa from interface to far region using normal direction
*******************************************************************************/
#ifdef DEBUG
  boil::oout<<"cipcsl2_curv_interface_ext: start\n";
#endif

  const real dtau = 0.25*dxmin;
  real flux,isgn;
  int mmax=8;  

  /*------------+
  |  set_wflag  |
  +------------*/
  set_wflag();

  /*-----------------------------------------------+
  |  initial value (necessary for high curvature)  |
  +-----------------------------------------------*/
  for(int flag = 1; flag <=6; flag++) {
    for_ijk(i,j,k){
      if (wflag[i][j][k]==flag) {
        real sumk = 0.0;
        int  sumc = 0;
        if (wflag[i-1][j][k]>-1000 && wflag[i-1][j][k]<=flag-1) {
          sumk += kappa[i-1][j][k];
          sumc ++;
        }
        if (wflag[i+1][j][k]>-1000 && wflag[i+1][j][k]<=flag-1) { 
          sumk += kappa[i+1][j][k];
          sumc ++;
        }
        if (wflag[i][j-1][k]>-1000 && wflag[i][j-1][k]<=flag-1) { 
          sumk += kappa[i][j-1][k];
          sumc ++;
        }
        if (wflag[i][j+1][k]>-1000 && wflag[i][j+1][k]<=flag-1) {
          sumk += kappa[i][j+1][k];
          sumc ++;
        }
        if (wflag[i][j][k-1]>-1000 && wflag[i][j][k-1]<=flag-1) {
          sumk += kappa[i][j][k-1];
          sumc ++;
        }
        if (wflag[i][j][k+1]>-1000 && wflag[i][j][k+1]<=flag-1) {
          sumk += kappa[i][j][k+1];
          sumc ++;
        }
        if (sumc == 0) {
          boil::aout<<"cipcsl2_curv_interface_ext: Error!!!\n";
          boil::aout<<"sumc=0 (increasing) flag= "<<flag<<"\n";
          boil::aout<<"i= "<<i<<" j= "<<j<<" k= "<<k<<"\n";
          exit(0);
        }
        fn[i][j][k] = sumk / real(sumc);
      }
    }

    /* update */
    for_ijk(i,j,k){
      if (wflag[i][j][k]==flag) {
        kappa[i][j][k]=fn[i][j][k];
      }
    }
    insert_bc_kappa(kappa);    // boundary condition
    kappa.exchange_all();
  }
  for(int flag = -1; flag >=-6; flag--) {
    for_ijk(i,j,k){
      if (wflag[i][j][k]==flag) {
        real sumk = 0.0;
        int  sumc = 0;
        if (wflag[i-1][j][k]>=flag+1) { 
          sumk += kappa[i-1][j][k];
          sumc ++;
        }
        if (wflag[i+1][j][k]>=flag+1) {
          sumk += kappa[i+1][j][k];
          sumc ++;
        }
        if (wflag[i][j-1][k]>=flag+1) {
          sumk += kappa[i][j-1][k];
          sumc ++;
        }
        if (wflag[i][j+1][k]>=flag+1) {
          sumk += kappa[i][j+1][k];
          sumc ++;
        }
        if (wflag[i][j][k-1]>=flag+1) {
          sumk += kappa[i][j][k-1];
          sumc ++;
        }
        if (wflag[i][j][k+1]>=flag+1) {
          sumk += kappa[i][j][k+1];
          sumc ++;
        }
        if (sumc == 0) {
          boil::aout<<"cipcsl2_curv_interface_ext: Error!!!\n";
          boil::aout<<"sumc=0 (decreasing)\n";
          exit(0);
        }
        fn[i][j][k] = sumk / real(sumc);
      }
    }

    /* update */
    for_ijk(i,j,k){
      if (wflag[i][j][k]==flag) {
        kappa[i][j][k]=fn[i][j][k];
      }
    }
    insert_bc_kappa(kappa);    // boundary condition
    kappa.exchange_all();
  }

  /*-----------------+
  |  solve equation  |
  +-----------------*/
  for(int mstep=1; mstep<=mmax; mstep++){
    /* advection term */
    for_ijk(i,j,k) {
      if(abs(wflag[i][j][k])>=1 && abs(wflag[i][j][k])<=6){
 
        isgn=copysign(1.0,dist[i][j][k]);
        fn[i][j][k]=0.0;
        /* i-direction */
        if (abs(wflag[i-1][j][k])>=7) {
          if (isgn * nx[i][j][k] > 0.0) {
            fn[i][j][k] -= isgn * (-nx[i][j][k])   // reverse 
                         * (kappa[i+1][j][k]-kappa[i][j][k])/dxe(i);
          } else {
            fn[i][j][k] -= isgn * nx[i][j][k] 
                         * (kappa[i+1][j][k]-kappa[i][j][k])/dxe(i);
          }
        } else if (abs(wflag[i+1][j][k])>=7) {
          if (isgn * nx[i][j][k] >0.0 ) {
            fn[i][j][k] -= isgn * nx[i][j][k] 
                         * (kappa[i][j][k]-kappa[i-1][j][k])/dxw(i);
          } else {
            fn[i][j][k] -= isgn * (-nx[i][j][k])   // reverse
                         * (kappa[i][j][k]-kappa[i-1][j][k])/dxw(i);
          }
        } else {
          if (isgn*nx[i][j][k] >=0 ){
            fn[i][j][k] -= isgn * nx[i][j][k] 
                         * (kappa[i][j][k]-kappa[i-1][j][k])/dxw(i);
          } else {
            fn[i][j][k] -= isgn * nx[i][j][k] 
                         * (kappa[i+1][j][k]-kappa[i][j][k])/dxe(i);
          }
        }

        /* j-direction */
        if (abs(wflag[i][j-1][k])>=7) {
          if (isgn * ny[i][j][k] > 0.0) {
            fn[i][j][k] -= isgn * (-ny[i][j][k])   // reverse
                         * (kappa[i][j+1][k]-kappa[i][j][k])/dyn(j);
          } else {
            fn[i][j][k] -= isgn * ny[i][j][k] 
                         * (kappa[i][j+1][k]-kappa[i][j][k])/dyn(j);
          }
        } else if (abs(wflag[i][j+1][k])>=7) {
          if (isgn * ny[i][j][k] > 0.0) {
            fn[i][j][k] -= isgn * ny[i][j][k] 
                         * (kappa[i][j][k]-kappa[i][j-1][k])/dys(j);
          } else {
            fn[i][j][k] -= isgn * (-ny[i][j][k])   // reverse 
                         * (kappa[i][j][k]-kappa[i][j-1][k])/dys(j);
          }
        } else {
          if (isgn*ny[i][j][k] >=0 ){
            fn[i][j][k] -= isgn * ny[i][j][k] 
                         * (kappa[i][j][k]-kappa[i][j-1][k])/dys(j);
          } else {
            fn[i][j][k] -= isgn * ny[i][j][k] 
                         * (kappa[i][j+1][k]-kappa[i][j][k])/dyn(j);
          }
        }

        /* k-direction */
        if (abs(wflag[i][j][k-1])>=7) {
          if (isgn * nz[i][j][k] > 0.0 ) {
            fn[i][j][k] -= isgn * (-nz[i][j][k])  // reverse 
                         * (kappa[i][j][k+1]-kappa[i][j][k])/dzt(k);
          } else {
            fn[i][j][k] -= isgn * nz[i][j][k] 
                         * (kappa[i][j][k+1]-kappa[i][j][k])/dzt(k);
          }
        } else if (abs(wflag[i][j][k+1])>=7) {
          if (isgn * nz[i][j][k] > 0.0) {
            fn[i][j][k] -= isgn * nz[i][j][k] 
                         * (kappa[i][j][k]-kappa[i][j][k-1])/dzb(k);
          } else {
            fn[i][j][k] -= isgn * (-nz[i][j][k])   // reverse 
                         * (kappa[i][j][k]-kappa[i][j][k-1])/dzb(k);
          }
        } else {
          if (isgn*nz[i][j][k] >=0 ){
            fn[i][j][k] -= isgn * nz[i][j][k] 
                         * (kappa[i][j][k]-kappa[i][j][k-1])/dzb(k);
          } else {
            fn[i][j][k] -= isgn * nz[i][j][k] 
                         * (kappa[i][j][k+1]-kappa[i][j][k])/dzt(k);
          }
        }
      }
    }

    /* update */
    for_ijk(i,j,k) {
      if(abs(wflag[i][j][k])>=1 && abs(wflag[i][j][k])<=6){
        kappa[i][j][k] += dtau*fn[i][j][k];
      }
    }

    //bnd_sym_kappa();       // must not call here
    insert_bc_kappa(kappa);  // boundary condition
    kappa.exchange_all();
  }

#ifdef DEBUG
  boil::oout<<"cipcsl2_curv_interface_ext:end\n";
#endif

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_curv_interface_ext.cpp,v 1.1 2015/05/06 08:00:09 sato Exp $'/
+-----------------------------------------------------------------------------*/
