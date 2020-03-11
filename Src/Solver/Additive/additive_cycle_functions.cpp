#include "additive.h"

//#define DEBUG

#ifdef DEBUG 
const char * sname = "solver"; 
#else
const char * sname = NULL;
#endif

/******************************************************************************/
bool AC::init_cycles(real & res_0, int * ncyc) {
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
  res_0 = residual(*L[0]);
  L[0]->fold = L[0]->phi; /* store "good" solution */
  if((res_0 < targ_res_val && min_cyc==0) || res_0==0.0) {
    boil::oout << "Converged in 0 cycles!" << boil::endl;
    if(ncyc) * ncyc = 0;
    return true;
  }
  
  boil::oout << "Initial res = " << res_0 << boil::endl;

  return false;
}

/******************************************************************************/
int AC::converged(const ResRat & factor, const int & cycle, 
                  const real & res_0, const real & res0,
                  int * ncyc) {
  real res1 = residual(*L[0]);

  boil::oout << "Cycle " << cycle << "; res = " << res1 << boil::endl;

  if( (res1/res_0 < factor || res1 < targ_res_val) && cycle >= min_cyc ) {
    if(cycle > 1)
      boil::oout << "Converged in " << cycle << " cycles!" << boil::endl;
    else
      boil::oout << "Converged in " << cycle << " cycle!" << boil::endl;

    return 0;
  }

  //if(stop_if_div)
  if(stop_if_div && cycle >= min_cyc)
    if(res1 >= res0) {
      L[0]->phi = L[0]->fold; /* restore last "good" solution */
      boil::oout << "Failed to conv. " << cycle << " cycles!" << boil::endl;
      return 1;
    }

  return 2;
}

/******************************************************************************/
void AC::call_smoother(const int l, const MaxIter & mi, 
                       const ResRat & res_rat, const ResTol & res_tol) {

  L[l] -> solver->solve( L[l]->A, L[l]->phi, L[l]->fnew,
                         mi, sname, res_rat, res_tol );

  return;
}

/******************************************************************************/
void AC::call_solver(const int l, const MaxIter & mi, 
                     const ResRat & res_rat, const ResTol & res_tol) {

  solver->solve( L[l]->A, L[l]->phi, L[l]->fnew,
                 mi, sname, res_rat, res_tol );

  return;
}
