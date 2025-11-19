#include "phasechangevof.h"

real PhaseChangeVOF::lambda(const int i, const int j, const int k,
                            const Scalar * diff_eddy) const {

  real lam;
  if(dom->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    if(clr[i][j][k]<clrsurf) {
      lam=lambdav;
      if(diff_eddy) lam += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
    } else {
      lam=lambdal;
      if(diff_eddy) lam += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
    }
  }

  return lam;
}
