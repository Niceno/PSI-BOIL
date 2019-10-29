#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::forward_cartesian(Scalar & scp) {
/***************************************************************************//**
 \brief Solve the forward Cartesian problem, i.e. calculate phi(alp,n)
    prerequisites: nx, nalpha
    output: tentative phi = scp
*******************************************************************************/

  for_vijk(scp,i,j,k) {
    real nnx  = -nx[i][j][k];
    real nny  = -ny[i][j][k];
    real nnz  = -nz[i][j][k];
    real nsum = fabs(nnx)+fabs(nny)+fabs(nnz);
    real n2sum = nnx*nnx+nny*nny+nnz*nnz;
    real nalp = nalpha[i][j][k];
    if       (!boil::realistic(nalp)) {
      if(nalp>0.0) {
        scp[i][j][k] = 1.0;
      } else {
        scp[i][j][k] = 0.0;
      }
    } else if(dom->ibody().off(i,j,k)|| n2sum <0.5) {
      scp[i][j][k] = -1.0;
    } else {
      /* normalization */
      nalp /= nsum;
      nnx  /= nsum;
      nny  /= nsum;
      nnz  /= nsum;
      scp[i][j][k] = calc_v(nalp,fabs(nnx),fabs(nny),fabs(nnz));
    }
  }

  scp.bnd_update();
  scp.exchange_all();

  return;
}


