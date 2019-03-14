#include "vof.h"

/******************************************************************************/
/* calculate liquid velocity */
/******************************************************************************/
void VOF::cal_liq_vel() {
  
  for_ijk(i,j,k) {

    /* calculate normal vector  
     * n points to the liquid
     * m is the real space normal vector  */
    real mmx = -mx[i][j][k];
    real mmy = -my[i][j][k];
    real mmz = -mz[i][j][k];

    /* cell centre velocity */
    real uxc = 0.5 * ((*u)[Comp::u()][i][j][k] + (*u)[Comp::u()][i+1][j][k]);
    real uyc = 0.5 * ((*u)[Comp::v()][i][j][k] + (*u)[Comp::v()][i][j+1][k]);
    real uzc = 0.5 * ((*u)[Comp::w()][i][j][k] + (*u)[Comp::w()][i][j][k+1]);

    /* normal velocity amplitude */
    real normvel = uxc*mmx+uyc*mmy+uzc*mmz;

    /* tangential velocity vector */
    uxc -= normvel*mmx;
    uyc -= normvel*mmy;
    uzc -= normvel*mmz;

    utx[i][j][k] = uxc;
    uty[i][j][k] = uyc;
    utz[i][j][k] = uzc;

    /* tangential velocity amplitude */
    real tangvel = sqrt(uxc*uxc+uyc*uyc+uzc*uzc);
 
    unliq[i][j][k] = normvel;
    utliq[i][j][k] = tangvel;
  }

  unliq.exchange();
  utliq.exchange();
  utx.exchange();
  uty.exchange();
  utz.exchange();

  //boil::plot->plot(unliq,utliq,utx,uty,utz, "un-ut-utx-uty-utz", time->current_step());

  return;
}
