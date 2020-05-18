#include "phasechange4.h"

real PhaseChange4::lambda(const int i, const int j, const int k,
                          const Scalar * diff_eddy) const {

  real lam;
  if(dom->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    if(!topo->above_interface(i,j,k)) {
      lam=lambdav;
      if(diff_eddy) lam += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
    } else {
      lam=lambdal;
      if(diff_eddy) lam += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
    }
  }

  return lam;
}

real PhaseChange4::lambda_inv(const int i, const int j, const int k,
                              const Scalar * diff_eddy) const {

  real lam;
  if(dom->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    /* note the inversion! */
    if(!topo->above_interface(i,j,k)) {
      lam=lambdal;
      if(diff_eddy) lam += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
    } else {
      lam=lambdav;
      if(diff_eddy) lam += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
    }
  }

  return lam;
}
