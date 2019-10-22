#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::backward_axisymmetric(const Scalar & vf, Scalar & alp) {
/***************************************************************************//**
 \brief Solve the backward axisymmetric problem, i.e. calculate alp(phi,n,x)
    prerequisites: nx
    output: alp = nalpha (= denormalized)
*******************************************************************************/

  for_vijk(alp,i,j,k) {
    /* n points to the liquid */
    real nnx  = -nx[i][j][k];
    real nnz  = -nz[i][j][k];
    real nsum = fabs(nnx)+fabs(nnz);
    real v = vf[i][j][k];
    real deltax = vf.dxc(i);
    if(  dom->ibody().off(i,j,k)|| nsum <boil::pico
       ||deltax==0.0) {
      alp[i][j][k] = ((v>=phisurf)-(v<phisurf))*boil::unreal;
    } else {
      /* normalization */
      nnx  /= nsum;
      real eta0 = vf.xn(i)/deltax;
      alp[i][j][k] = nsum*calc_alpha_axisymmetric(nnx,v,eta0);
    }
  }

  alp.bnd_update();
  alp.exchange_all();

  return;
}


