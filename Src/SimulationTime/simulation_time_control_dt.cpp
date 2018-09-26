#include "simulation_time.h"

/******************************************************************************/
void Times::control_dt(const real current_cfl,
                       const real target_cfl,
                       const real dt_max) {

  // coefficient for increase dt: coef_inc
  // coefficient for decrease dt: coef_dec

  /* increase nhist */
  increase_nhist();

  /* shift history */
  int nhist = get_nhist();
  if(nhist==CFL_HIST_SIZE){
    for(int ihist=1; ihist<=nhist-1; ihist++){
      set_cfl_hist(ihist-1, get_cfl_hist(ihist));
    }
  }

  /* store current cfl */
  set_cfl_hist(nhist-1, current_cfl);

  real d_t = dt();
  if (current_cfl>=target_cfl) {
    /* decrease dt */
    real dt_new = std::min(dt_max, target_cfl/current_cfl*d_t);
    dt_new = std::max(coef_dec*d_t, dt_new);
    set_dt(dt_new);
    boil::oout << "@control_dt; decrease dt: dt = " << dt_new << boil::endl;
  } else if (d_t>dt_max) {
    real dt_new  = std::max(coef_dec*d_t , dt_max);
    set_dt(dt_new);
    boil::oout << "@control_dt; decrease dt: dt = " << dt_new << boil::endl;
  } else if (nhist==CFL_HIST_SIZE) {
    /* check history */
    int iexceed=0;
    for(int ihist=0; ihist<=nhist-1; ihist++){
      if (get_cfl_hist(ihist)>=target_cfl) iexceed=1;
    }
    /* increase dt */
    if (iexceed==0 && current_step()>=istep_inc+8) {
      real dt_new = std::min(coef_inc*d_t, 
                             target_cfl/(current_cfl+boil::pico)*d_t);
      dt_new = std::min(dt_new,dt_max);
      if (dt_new>d_t) {
        set_dt(dt_new);
        istep_inc = current_step();
        boil::oout << "@control_dt; increase dt: dt = " << dt_new << boil::endl;
      }
    }
  }

  if (nhist==CFL_HIST_SIZE) decrease_nhist(); 

}

/*-----------------------------------------------------------------------------+
 '$Id: simulation_time_control_dt.cpp,v 1.6 2015/06/29 18:17:56 sato Exp $'/
+-----------------------------------------------------------------------------*/
