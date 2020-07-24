#include "integralpc.h"

/******************************************************************************/
real IntegralPC::regula_falsi_kernel(real & tm, real & tp, 
                                     real & vpcm, real & vpcp,
                                     real & vpc_tpr_m, real & vpc_tpr_p,
                                     real & vpc_cav_m, real & vpc_cav_p,
                                     const real pinf) {
/***************************************************************************//**
*  \brief root finding algorithm kernel. 
*******************************************************************************/

  /* regula falsi */
  real t_next = (tm*vpcp-tp*vpcm)/(vpcp-vpcm);

  real vpc_tpr_next,vpc_cav_next;
  real vpc_next = evaluate_point(t_next,pinf,
                                 vpc_tpr_next, vpc_cav_next);

  if(vpc_next*vpcm>0.) {
    tm = t_next;
    vpc_tpr_m = vpc_tpr_next;
    vpc_cav_m = vpc_cav_next;
    vpcm = vpc_next;
  } else {
    tp = t_next;
    vpc_tpr_p = vpc_tpr_next;
    vpc_cav_p = vpc_cav_next;
    vpcp = vpc_next;
  }

  real eps_calc = std::min(fabs(vpcm/std::min(fabs(vpc_tpr_m),fabs(vpc_cav_m)))
                          ,fabs(vpcp/std::min(fabs(vpc_tpr_p),fabs(vpc_cav_p))));

  return eps_calc;
}
