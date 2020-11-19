#include "vof.h"

/******************************************************************************/
void VOF::smooth(Scalar & sca, const int itnum) {
/***************************************************************************//**
* \brief  Smooth/sharpen scalar (curvature).
*            input    : sca, itnum
*            output   : sca
*            temporary: stmp
*******************************************************************************/

  stmp = 0.;

  if (itnum>=1) {
    /*-----------------------------+
    |  diffusion, explicit Jakobi  |
    +-----------------------------*/
    real dtau=0.125;
    /* iterate */
    for(int it=0; it<itnum; it++) {
      for_ijk(i,j,k) {
        if(dom->ibody().off(i,j,k)) continue;
        if(!boil::realistic(sca[i][j][k])) continue;

        real coef_x_m = bflag_struct.ifull ? 1.0 : 0.0;
        real coef_x_p = bflag_struct.ifull ? 1.0 : 0.0;
        real coef_y_m = bflag_struct.jfull ? 1.0 : 0.0;
        real coef_y_p = bflag_struct.jfull ? 1.0 : 0.0;
        real coef_z_m = bflag_struct.kfull ? 1.0 : 0.0;
        real coef_z_p = bflag_struct.kfull ? 1.0 : 0.0;

        if(!boil::realistic(sca[i-1][j][k]))
          coef_x_m = 0.0;
        if(!boil::realistic(sca[i+1][j][k]))
          coef_x_p = 0.0;
        if(!boil::realistic(sca[i][j-1][k]))
          coef_y_m = 0.0;
        if(!boil::realistic(sca[i][j+1][k]))
          coef_y_p = 0.0;
        if(!boil::realistic(sca[i][j][k-1]))
          coef_z_m = 0.0;
        if(!boil::realistic(sca[i][j][k+1]))
          coef_z_p = 0.0;

        real coef_c = coef_x_m + coef_x_p
                    + coef_y_m + coef_y_p
                    + coef_z_m + coef_z_p;
        if(coef_c<0.5) continue;

        stmp[i][j][k] = coef_x_m*sca[i-1][j][k] + coef_x_p*sca[i+1][j][k]
                      + coef_y_m*sca[i][j-1][k] + coef_y_p*sca[i][j+1][k]
                      + coef_z_m*sca[i][j][k-1] + coef_z_p*sca[i][j][k+1]
                      - coef_c*sca[i][j][k];
      } /* ijk */

      /* update */
      for_ijk(i,j,k) {
        sca[i][j][k]+=dtau*stmp[i][j][k];
      }
      sca.bnd_update();
      sca.exchange();
    } /* iter */
  } /* itnum >= 1 */

#if 0
  boil::plot->plot(sca, "sca", time->current_step());
  exit(0);
#endif

  return;
}
