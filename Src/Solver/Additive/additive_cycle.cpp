#include "additive.h"

/******************************************************************************/
bool AC::cycle(const Cycle & init, const Cycle & loop, const ResTol & toler,
               const ResRat & factor, const std::array<MaxIter,3> & mv,
               const std::array<MaxIter,3> & ms, int * ncyc) {
/*------------------+
|  make few cycles  |
+------------------*/

  boil::timer.start("cycle");

  /*--------------------+
  |  initialize cycles  |
  +--------------------*/
  real res_0, reslinf_0;
  if(init_cycles(toler,res_0,reslinf_0,ncyc)) {
    boil::oout<<"res_0 "<<res_0<<boil::endl;
    boil::timer.stop("cycle");
    return true;
  }

  /*-----------------------------------------------------+
  |    very important: initialize all levels to zero!    |
  |  without it, algorithm has difficulties to converge  | 
  |    after a certain, (large) number of time steps     |
  |           (different for full multigrid)             |
  +-----------------------------------------------------*/
  for(int l=1; l<=nlevels-1; l++) {
    L[l]->phi  = 0.0;
    L[l]->fnew = 0.0;
  }

  /* convergence metrics */
  real reslinf_m0, res_m0, reslinf_m1, res_m1;

  if(init!=Cycle::none()) {
    //L[0]->phi = 0.0;
    full_cycle(0,init,mv,ms,0);

    /* evaluate convergence */
    res_m1 = residual(*L[0],&reslinf_m1);
    converged(toler,factor,0,res_m0,reslinf_m0,
                             res_0,reslinf_0,res_m1,reslinf_m1,ncyc);
  } else {
  
    /* evaluate convergence */
    res_m1 = residual(*L[0],&reslinf_m1);
  }

  /*--------------------------+
  |  loop through the cycles  |
  +--------------------------*/

  for(int c=1; c<=max_cyc; c++) {
 
    /* recursive kernel */
    if       (loop==Cycle::V()) {
      v_cycle(0,mv,ms,c);
    } else if(loop==Cycle::F()) {
      f_cycle(0,mv,ms,c);
    } else if(loop==Cycle::W()) {
      w_cycle(0,mv,ms,c);
    } else if(loop==Cycle::flex()) {
      flex_cycle(0,mv,ms,c);
    /* just solve */
    } else {
      call_solver(0,MaxIter(40),ResRat(0.001),ResTol(boil::femto),MaxIter(-1),c);
      res_m0 = residual(*L[0]);
      boil::oout << "res = " << res_m0 << boil::endl;
      break;
    }

    /* evaluate convergence */
    int con = converged(toler,factor,c,res_m0,reslinf_m0,
                                       res_0,reslinf_0,res_m1,reslinf_m1,ncyc);
    res_m1 = res_m0;
    reslinf_m1 = reslinf_m1;

    switch(con) {
      case 0 :
        boil::timer.stop("cycle");
        if(ncyc) * ncyc = c;
        return true;

      case 1 :
        boil::timer.stop("cycle");
        if(ncyc) * ncyc = c;
        return false;

      default :
        L[0]->fold = L[0]->phi; /* store the "good" solution */
    }
  }

  boil::timer.stop("cycle");
  if(ncyc) * ncyc = max_cyc;
  return false;
}
