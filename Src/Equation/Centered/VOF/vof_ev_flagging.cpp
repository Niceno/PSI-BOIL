#include "vof.h"

/******************************************************************************/
void VOF::ev_flagging(const Scalar & scp, const ScalarInt & iflag,
                      ScalarInt & otpflag, const Sign & sig) {
/***************************************************************************//**
*  \brief Flag the interface for the Poisson problem.
*  otpflag = -1 : color function <  clrsurf & at interface
*  otpflag = +1 : color function >= clrsurf & at interface 
*  otpflag = -2 : color function <  clrsurf & next to -1||1 cell
*  otpflag = +2 : color function >= clrsurf & next to -1||1 cell
*  otpflag = -3 : color function <  clrsurf & bulk
*  otpflag = +4 : color function >= clrsurf & bulk
*  otpflag = -1000: solid cell
*
*  sig == neg -> inversion of above
*******************************************************************************/
 
  otpflag = -1000;

  for_aijk(i,j,k) {
    if(  abs(iflag[i][j][k])==1 /* cells at the interface */
       ||iflag[i][j][k]==-3*sig) /* bulk gas cell */
      otpflag[i][j][k] = iflag[i][j][k]*sig;
    if(iflag[i][j][k]==3*sig) /* bulk liquid cell */
      otpflag[i][j][k] = 4;
  }

  for_ijk(i,j,k) {
    if(abs(otpflag[i][j][k])!=1) { /* cells next to the interface */
      if(  abs(otpflag[i-1][j][k])==1
         ||abs(otpflag[i+1][j][k])==1
         ||abs(otpflag[i][j-1][k])==1
         ||abs(otpflag[i][j+1][k])==1
         ||abs(otpflag[i][j][k-1])==1
         ||abs(otpflag[i][j][k+1])==1)
        otpflag[i][j][k] = 2*sig
                            *((scp[i][j][k]>phisurf)-(scp[i][j][k]<phisurf));
    }
  }

  for_aijk(i,j,k) {
    if(dom->ibody().off(i,j,k)) /* ibody */
      otpflag[i][j][k] = -1000;
  }

  for(int b=0; b<scp.bc().count(); b++) { /* inside wall */
    if(scp.bc().type(b) == BndType::wall()) {
      Dir d = scp.bc().direction(b);
      if(d != Dir::undefined()) {
        for_vijk(scp.bc().at(b), i,j,k) {
          otpflag[i][j][k]= -1000;
        }
      }
    }
  }

  otpflag.bnd_update_nowall();
  otpflag.exchange();

  return;
}
