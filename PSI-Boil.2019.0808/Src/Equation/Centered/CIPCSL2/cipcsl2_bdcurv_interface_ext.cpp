#include "cipcsl2.h"
#include "../../../Parallel/Out/out.h"
#include <iomanip>
using namespace std;

/******************************************************************************/
void CIPCSL2::bdcurv_interface_ext() {
/***************************************************************************//**
*  \brief extrapolate kappa to normal direction in wall adjacent layer
*  crude code: need to extrapolate more accurately
*******************************************************************************/
#ifdef DEBUG
  boil::oout<<"cipcsl2_ext_kappa:start\n";
#endif

#if 0
  cout<<"bdcurv_int_ext: "<<kappa[4][17][1]<<" "<<kappa[3][17][1]<<"\n";
#endif
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

#if 0
  boil::plot->plot(clr, wflag, "clr-wflag", time->current_step());
  boil::plot->plot(clr, kappa, "clr-kappa", time->current_step());
  exit(0);
#endif

#ifdef DEBUG
  boil::oout<<"cipcsl2_bdcurv_interface_ext:end\n";
#endif

  return;
}
