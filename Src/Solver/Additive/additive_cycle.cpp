#include "additive.h"

/******************************************************************************/
bool AC::cycle(const Cycle & init, const Cycle & loop,
               const ResRat & factor, int * ncyc) {
/*------------------+
|  make few cycles  |
+------------------*/

  boil::timer.start("cycle");

  /*--------------------+
  |  initialize cycles  |
  +--------------------*/
  real res_0;
  if(init_cycles(res_0,ncyc)) {
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
  for(int l=1; l<nlevels-1; l++) {
    L[l]->phi  = 0.0;
    L[l]->fnew = 0.0;
  }
  if(init!=Cycle::none()) {
    L[0]->phi = 0.0;
    L[nlevels-1]->fnew = 0.0;
    full_cycle(0,init);

    /* evaluate convergence */
    real res0 = residual(*L[0]);
    int con = converged(factor,0,res_0,res0,ncyc);
  } else {
    L[nlevels-1]->phi  = 0.0;
    L[nlevels-1]->fnew = 0.0;
  }

  /*--------------------------+
  |  loop through the cycles  |
  +--------------------------*/
  for(int c=1; c<=max_cyc; c++) {
 
    real res0 = residual(*L[0]);

    /* recursive kernel */
    if       (loop==Cycle::V()) {
      v_cycle(0);
    } else if(loop==Cycle::F1()) {
      f_cycle(0,false);
    } else if(loop==Cycle::F2()) {
      f_cycle(0,true);
    } else if(loop==Cycle::W1()) {
      w_cycle(0,false);
    } else if(loop==Cycle::W2()) {
      w_cycle(0,true);
    /* just solve */
    } else {
      call_solver(0,MaxIter(40),ResRat(0.001),ResTol(boil::femto));
      real res1 = residual(*L[0]);
      boil::oout << "res = " << res1 << boil::endl;
      break;
    }

    /* evaluate convergence */
    int con = converged(factor,c,res_0,res0,ncyc);
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
