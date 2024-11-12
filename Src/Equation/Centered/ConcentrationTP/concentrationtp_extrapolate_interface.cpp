#include "concentrationtp.h"

/******************************************************************************/
void ConcentrationTP::extrapolate_interface() {
/***************************************************************************//**
*  \brief Extrapolate eps in the cellc which newly appear in vapor phase
*******************************************************************************/

#if 0
  /* visualize eps */
  boil::plot->plot(clr,phi,"before-c-eps", time->current_step());
  //exit(0);
#endif

  stmp=phi;

  for_ijk(i,j,k) {
    real col_new = vfval(i,j,k);
    real col_old = vfvalold(i,j,k);

    // if color function is the volume fraction of liquid
    if(matter_sig==Sign::neg()) col_new = 1.0-col_new;
    if(matter_sig==Sign::neg()) col_old = 1.0-col_old;

    // detect the cells, which was full-of-liquid in the last step
    // but is including vapor in the current step
    if( col_old < col_crit
     && col_new > col_crit) {
      //std::cout<<"gotcha! "<<i<<" "<<j<<" "<<k<<"\n";

      // extrapolate from surrouding cells
      // weight
      real wim=0.0, wip=0.0, wjm=0.0, wjp=0.0, wkm=0.0, wkp=0.0;

      // i-
      real col_old_im = vfvalold(i-1,j,k);
      if(matter_sig==Sign::neg()) col_old_im = 1.0-col_old_im;
      // i- included vapor in the previous step
      if( col_old_im > col_crit ) { wim=1.0; }

      // i+
      real col_old_ip = vfvalold(i+1,j,k);
      if(matter_sig==Sign::neg()) col_old_ip = 1.0-col_old_ip;
      // i+ included vapor in the previous step
      if( col_old_ip > col_crit ) { wip=1.0; }

      // j-
      real col_old_jm = vfvalold(i,j-1,k);
      if(matter_sig==Sign::neg()) col_old_jm = 1.0-col_old_jm;
      // j- included vapor in the previous step
      if( col_old_jm > col_crit ) { wjm=1.0; }

      // j+
      real col_old_jp = vfvalold(i,j+1,k);
      if(matter_sig==Sign::neg()) col_old_jp = 1.0-col_old_jp;
      // j+ included vapor in the previous step
      if( col_old_jp > col_crit ) { wjp=1.0; }

      // k-
      real col_old_km = vfvalold(i,j,k-1);
      if(matter_sig==Sign::neg()) col_old_km = 1.0-col_old_km;
      // k- included vapor in the previous step
      if( col_old_km > col_crit ) { wkm=1.0; }

      // k+
      real col_old_kp = vfvalold(i,j,k+1);
      if(matter_sig==Sign::neg()) col_old_kp = 1.0-col_old_kp;
      // j+ included vapor in the previous step
      if( col_old_kp > col_crit ) { wkp=1.0; }

      real w_sum = wim+wip+wjm+wjp+wkm+wkp;
      if (w_sum > 0) {
      stmp[i][j][k] = (wim*phi[i-1][j][k] + wip*phi[i+1][j][k]
                    +  wjm*phi[i][j-1][k] + wjp*phi[i][j+1][k]
                    +  wkm*phi[i][j][k-1] + wkp*phi[i][j][k+1])/w_sum;
      }
    }
  }

  stmp.bnd_update();
  stmp.exchange();
  phi = stmp;

  #if 0
    /* visualize eps */
    boil::plot->plot(clr,phi,"after-c-eps", time->current_step());
    exit(0);
  #endif

  return;
}
