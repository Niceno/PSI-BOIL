#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::forward_axisymmetric(const Scalar & color,
                                     Scalar & axip, Scalar & Kp) {
/***************************************************************************//**
 \brief Solve the forward axisymmetric problem, i.e. calculate phi(alp,n,x)
    prerequisites: nx, nalpha, clr
    output: tentative phi = axip and the correction coefficient K = phi/c
*******************************************************************************/

  //boil::plot->plot(color,nx,ny,nz,nalpha,"clr-nx-ny-nz-nalp", time->current_step());

  for_vijk(axip,i,j,k) {
    /* n points to the liquid */
    real nnx  = -nx[i][j][k];
    real nnz  = -nz[i][j][k];
    real nsum = fabs(nnx)+fabs(nnz);
    real nalp = nalpha[i][j][k];
    real deltax = axip.dxc(i);
    if(  dom->ibody().off(i,j,k)|| nsum <boil::pico
       ||!boil::realistic(nalp) ||deltax==0.0) {
      axip[i][j][k] = color[i][j][k];
      Kp[i][j][k] = 1.0;
    } else {
      /* normalization */
      nalp /= nsum;
      nnx  /= nsum;
      real eta0 = axip.xn(i)/deltax;
      axip[i][j][k] = calc_v_axisymmetric(nnx,nalp,eta0, Kp[i][j][k]);
      //boil::oout<<i<<" "<<j<<" "<<k<<" | "<<Kp[i][j][k]<<boil::endl;
    }
  }

  axip.bnd_update();
  axip.exchange_all();
  Kp.bnd_update();
  Kp.exchange_all();

  return;
}
