#include "commonheattransfer.h"

real CommonHeatTransfer::lambda(const int i, const int j, const int k,
                                const Scalar * diff_eddy) const {

  real lam;
  if(topo->domain()->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    if(!topo->above_interface(i,j,k)) {
      lam=lambdav(i,j,k,diff_eddy);
    } else {
      lam=lambdal(i,j,k,diff_eddy);
    }
  }

  return lam;
}

real CommonHeatTransfer::lambda_inv(const int i, const int j, const int k,
                                    const Scalar * diff_eddy) const {

  real lam;
  if(topo->domain()->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    /* note the inversion! */
    if(!topo->above_interface(i,j,k)) {
      lam=lambdal(i,j,k,diff_eddy);
    } else {
      lam=lambdav(i,j,k,diff_eddy);
    }
  }

  return lam;
}

real CommonHeatTransfer::lambda_old(const int i, const int j, const int k,
                                    const Scalar * diff_eddy) const {

  real lam;
  if(topo->domain()->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    if(!topo->above_interface_old(i,j,k)) {
      lam=lambdav(i,j,k,diff_eddy);
    } else {
      lam=lambdal(i,j,k,diff_eddy);
    }
  }

  return lam;
}

real CommonHeatTransfer::lambda_inv_old(const int i, const int j, const int k,
                                        const Scalar * diff_eddy) const {

  real lam;
  if(topo->domain()->ibody().off(i,j,k)) {
    lam = solid()->lambda(i,j,k);
  } else {
    /* note the inversion! */
    if(!topo->above_interface_old(i,j,k)) {
      lam=lambdal(i,j,k,diff_eddy);
    } else {
      lam=lambdav(i,j,k,diff_eddy);
    }
  }

  return lam;
}
