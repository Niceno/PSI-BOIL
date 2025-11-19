#include "vof.h"

/******************************************************************************/
void VOF::color_minmax() {
/***************************************************************************//**
*  \brief Detect maximum and minimum of color function.
*           results: cmin, cmax, ijk
*******************************************************************************/
  real cmin = +boil::exa;
  real cmax = -boil::exa;

  int i_cmin, j_cmin, k_cmin, i_cmax, j_cmax, k_cmax;
  i_cmin = j_cmin = k_cmin = i_cmax = j_cmax = k_cmax = -1;

  for_ijk(i,j,k){
    if ( dom->ibody().on(i,j,k) ) {
      if( cmin > phi[i][j][k] ) {
        cmin = phi[i][j][k];
        i_cmin = i;
        j_cmin = j;
        k_cmin = k; 
      }
      if( cmax < phi[i][j][k] ) {
        cmax = phi[i][j][k];
        i_cmax = i;
        j_cmax = j;
        k_cmax = k; 
      }
    }
  }
  real cminLocal = cmin;
  real cmaxLocal = cmax;

  boil::cart.min_real(&cmin);
  boil::cart.max_real(&cmax);

  set_minval(cmin);
  set_maxval(cmax);

  boil::oout<<"color_minmax: "<<time->current_time()<<" "
            <<cmin<<" "<<cmax<<"\n";

  if (cmin<-1.0 || cmax>2.0) {
    boil::oout<<"cipcsl2_color_minmax: Stop calculation because of";
    boil::oout<<" too small or too large color function.\n";
    if (cminLocal<-1.0 || cmaxLocal>2.0) { 
      boil::aout<<"proc= "<<boil::cart.iam()
        <<" min= "<<cminLocal<<" max= "<<cmaxLocal
        <<" i_min= "<<i_cmin-si()+1<<" j_min= "<<j_cmin-sj()+1
        <<" k_min= "<<k_cmin-sk()+1
        <<" i_max= "<<i_cmax-si()+1<<" j_max= "<<j_cmax-sj()+1
        <<" k_min= "<<k_cmax-sk()+1<<"\n";
    }
    boil::plot->plot((*u),phi,kappa, "uvw-phi-kappa", time->current_step());
    exit(0);
  }
}
