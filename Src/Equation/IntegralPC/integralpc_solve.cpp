#include "integralpc.h"
#include <iomanip>

/******************************************************************************/
real IntegralPC::solve(const ResTol & rt, const real pinf,
                       const Range<real> & tprr, const bool progress) {
/***************************************************************************//**
*  \brief Find new interfacial temperature.
*******************************************************************************/

  /* only boiling is considered atm */
  t_old = std::max(t_old,t_new);

  /* test old value */
  real vpc_tpr_old, vpc_cav_old;
  real vpc_old = evaluate_point(t_old,pinf,
                                vpc_tpr_old, vpc_cav_old);

  /* relative error */
  real eps_calc = fabs(vpc_old/vpc_cav_old);
  boil::oout<<"@@IPCbegin= "<<time.current_time();
  std::cout<< std::setprecision(16);
  boil::oout<<" "<<t_old<<" ";
  std::cout<< std::setprecision(8);
  boil::oout<<vpc_old<<" "<<vpc_cav_old<<" "<<vpc_tpr_old<<" "<<eps_calc<<boil::endl;

  /* find new value */
  t_new = t_old;
  if(eps_calc>=rt) {

    /* find bounds: dangerous crude code! */
    real mult = vpc_old > 0. ? -1. : +1.;
    real delta_t = 0.01;

    real vpc_tpr_new, vpc_cav_new, vpc_new;
    do {
      t_new = std::max( tprr.first(),
                        std::min(tprr.last(),t_new + mult * delta_t) );
      vpc_new = evaluate_point(t_new,pinf,
                               vpc_tpr_new, vpc_cav_new);
    } while(vpc_new*vpc_old>0.0);

    /* find root */
    int niter(0);
    eps_calc = std::min(fabs(vpc_old/vpc_cav_old),fabs(vpc_new/vpc_cav_new));
    boil::oout<<"@@IPCbrackets= "<<time.current_time()<<" "<<niter<<" ";
    std::cout<< std::setprecision(16);
    boil::oout<<t_old<<" "<<t_new<<" ";
    std::cout<< std::setprecision(8);
    boil::oout<<eps_calc
              <<" | "<<vpc_old<<" "<<vpc_cav_old<<" "<<vpc_tpr_old
              <<" | "<<vpc_new<<" "<<vpc_cav_new<<" "<<vpc_tpr_new<<boil::endl;

    while(eps_calc>=rt&&fabs(t_old-t_new)>errtol) {
      eps_calc = regula_falsi_kernel(t_old, t_new, vpc_old, vpc_new,
                                     vpc_tpr_old, vpc_tpr_new, 
                                     vpc_cav_old, vpc_cav_new,
                                     pinf);

      niter++;
      boil::oout<<"@@IPCvpc= "<<time.current_time()<<" "<<niter;
      std::cout<< std::setprecision(16);
      boil::oout<<" "<<t_old<<" "<<t_new<<" ";
      std::cout<< std::setprecision(8);
      boil::oout<<eps_calc
                <<" | "<<vpc_old<<" "<<vpc_cav_old<<" "<<vpc_tpr_old
                <<" | "<<vpc_new<<" "<<vpc_cav_new<<" "<<vpc_tpr_new<<boil::endl;
    }

    /* choose the better value */
    t_old = fabs(vpc_old/std::min(fabs(vpc_tpr_old),fabs(vpc_cav_old)))
          > fabs(vpc_new/std::min(fabs(vpc_tpr_new),fabs(vpc_cav_new)))
          ? t_new : t_old;
  }

  /* progress the actual time step? */
  if(progress) {
    /* restore fields */
    tpr = tpr_old;
    for_m(m)
      uvw(m) = uvw_old(m);

    /* standard phase-change */
    phase_change_evaluation(t_old);
    ns.vol_phase_change(&f);

    /* momentum */
    momentum_evaluation();
  }

  return t_old;
}
