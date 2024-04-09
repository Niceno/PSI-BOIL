#include "additive.h"

/******************************************************************************/
bool AC::vcycle(const ResRat & factor, int * ncyc) {
/*--------------------+
|  make few V-cycles  |
+--------------------*/

  boil::timer.start("vcycle");
	
  /*--------------------+ 
  |  update the r.h.s.  |
  +--------------------*/
  L[0]->update_rhs();

  /*------------------------------+ 
  |     update coarser levels     |
  + - - - - - - - - - - - - - - - +
  |   this creates overhead for   |
  |  systems which do not change  |
  +------------------------------*/
  for(int l=1; l<nlevels; l++) 
    coarsen_system(*L[l-1], *L[l]); // finer, coarser

  /*-------------------+
  |  initial residual  |
  +-------------------*/ 
  real res_0 = residual(*L[0]);
  L[0]->fold = L[0]->phi; /* store "good" solution */
  if((res_0 < targ_res_val && min_cyc==0) || res_0==0.0) {
    boil::oout << "Converged in 0 cycles!" << boil::endl;
    boil::timer.stop("vcycle");
    if(ncyc) * ncyc = 0;
    return true;
  }
  
  boil::oout << "Initial  res = " << res_0 << boil::endl;

  /*-----------------------------------------------------+
  |    very important: initialize all levels to zero!    |
  |  without it, algorithm has difficulties to converge  | 
  |    after a certain, (large) number of time steps     |
  +-----------------------------------------------------*/
  for(int l=1; l<nlevels; l++) {
    L[l]->phi  = 0.0;
    L[l]->fnew = 0.0;
  }

  /*--------------------------+
  |  loop through the cycles  |
  +--------------------------*/
  for(int cycle=1; cycle<=max_cyc; cycle++) {
 
    real res0 = residual(*L[0]);

    /* down */
    for(int l=0; l<nlevels-1; l++) {
      L[l] -> solver->solve( L[l]->A, L[l]->phi, L[l]->fnew, 
                             MaxIter(5), NULL, 
                             ResRat(0.1), ResTol(boil::nano) );

      residual(*L[l]);
      restriction(*L[l], *L[l+1]);
      //~L[l] -> print();
    }	    

    /* bottom and up */
    for(int l=nlevels-1; l>=0; l--) {
      if(l==nlevels-1) 
        L[l] -> solver->solve( L[l]->A, L[l]->phi, L[l]->fnew, 
                               MaxIter(40 * (l+1)), NULL, 
                               ResRat(0.001), ResTol(boil::femto) );
      else
        L[l] -> solver->solve( L[l]->A, L[l]->phi, L[l]->fnew, 
                               MaxIter(20 * (l+1)), NULL, 
                               ResRat(0.01), ResTol(boil::pico) );

      if(l > 0)
        interpolation(*L[l], *L[l-1]);
      //~L[l] -> print();
    }	    

    real res1 = residual(*L[0]);

    boil::oout << "Cycle " << cycle << "; res = " << res1 << boil::endl;
    
    if( (res1/res_0 < factor || res1 < targ_res_val) && cycle >= min_cyc ) {
      if(cycle > 1)
        boil::oout << "Converged in " << cycle << " cycles!" << boil::endl;
      else
        boil::oout << "Converged in " << cycle << " cycle!" << boil::endl;
//      L[0] -> print();
      boil::timer.stop("vcycle");
      if(ncyc) * ncyc = cycle;
      return true;
    }  

    //if(stop_if_div)
    if(stop_if_div && cycle >= min_cyc)
      if(res1 >= res0) {
        L[0]->phi = L[0]->fold; /* restore last "good" solution */
        boil::oout << "Failed to conv. " << cycle << " cycles!" << boil::endl;
//        L[0] -> print();
        boil::timer.stop("vcycle");
        if(ncyc) * ncyc = cycle;
        return false;
      }  
    L[0]->fold = L[0]->phi; /* store the "good" solution */
  }

  boil::timer.stop("vcycle");
  if(ncyc) * ncyc = max_cyc;
  return false;
}
