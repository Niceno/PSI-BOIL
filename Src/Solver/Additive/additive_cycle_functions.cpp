#include "additive.h"

//#define DEBUG

#ifdef DEBUG 
const char * sname = "solver"; 
#else
const char * sname = NULL;
#endif

/******************************************************************************/
bool AC::init_cycles(const ResTol & toler, real & res_0, real & reslinf_0,
                     int * ncyc) {
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
  res_0 = residual(*L[0],&reslinf_0);
  real & res_0_control = use_linf ? reslinf_0 : res_0;
  
  boil::oout << "Initial res = " << res_0 <<" ; "<<reslinf_0<< boil::endl;

  L[0]->fold = L[0]->phi; /* store "good" solution */
  if((res_0_control < real(toler) && min_cyc==0) || res_0==0.0) {
    boil::oout << "Converged in 0 cycles!" << boil::endl;
    if(ncyc) * ncyc = 0;
    return true;
  }

  return false;
}

/******************************************************************************/
int AC::converged(const ResTol & toler, const ResRat & factor,
                  const int & cycle, 
                  real & res1, real & reslinf1,
                  const real & res_beg, const real & reslinf_beg,
                  const real & res0, const real & reslinf0,
                  int * ncyc) {
  res1 = residual(*L[0],&reslinf1);

  boil::oout << "Cycle " << cycle << "; res = " << res1 <<" ; "
             << reslinf1<< boil::endl;

  real & res_control = use_linf ? reslinf1 : res1;
  const real & res_0_control = use_linf ? reslinf0 : res0;
  const real & res_beg_control = use_linf ? reslinf_beg : res_beg;

  if( (cycle >= min_cyc) && (   (res_control/res_beg_control < real(factor))
                             || (res_control < real(toler))  ) ) {
    if(cycle > 1)
      boil::oout << "Converged in " << cycle << " cycles!" << boil::endl;
    else
      boil::oout << "Converged in " << cycle << " cycle!" << boil::endl;

    return 0;
  }

  bool check_if_div=stop_if_div;
  if (priority_min_cyc) {
    if (cycle < min_cyc) {
      check_if_div=false;
    }
  }
  if (check_if_div) {
    if (res_control >= res_0_control) {
      L[0]->phi = L[0]->fold; /* restore last "good" solution */
      boil::oout << "Failed to conv. " << cycle << " cycles!" << boil::endl;
      return 1;
    }
  }
   
  return 2;
}

/******************************************************************************/
bool AC::call_smoother(const int l, const MinIter & mini, const MaxIter & mi, 
                       const ResRat & res_rat, const ResTol & res_tol,
                       const MaxIter & ms, const int gi) {

  return L[l]->solver->solve(L[l]->A, L[l]->phi, L[l]->fnew,
                             mini, mi, sname, res_rat, res_tol, 
                             L[0]->time->dti()*L[0]->scale, int(ms));
}

/******************************************************************************/
bool AC::call_solver(const int l, const MinIter & mini, const MaxIter & mi, 
                     const ResRat & res_rat, const ResTol & res_tol,
                     const MaxIter & ms, const int gi) {

  return solver->solve(L[l]->A, L[l]->phi, L[l]->fnew,
                       mini, mi, sname, res_rat, res_tol, 
                       L[0]->time->dti()*L[0]->scale, int(ms));
}
