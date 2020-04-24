#include "additive.h"

/***************************************************************************//**
*  Constructs the AC solver. As a first step, coarsend the Centered variable
*  for which it is meant to be used. As a second step, stores pointers to
*  variable's levels.  
*******************************************************************************/
AC::AC(Centered * cen, Linear * sol) {

  /* set variables which steer the cycle */
  max_cyc = 20;
  min_cyc =  0;
  targ_res_val = boil::nano;
  targ_res_rat = 0.01;
  stop_if_div  = true;
  mv_def = {MaxIter(20),MaxIter(20),MaxIter(20)};

  /* coarsen the variable */
  cen -> coarsen();	
  
  /* store the pointers to variable levels */
  L[0] = cen;
  for(int i=1;   ;i++) {
    if( L[i-1]->coarser() != NULL ) {
      L[i] = L[i-1]->coarser(); 
    } else {
      nlevels = i;
      break;
    }
  }

  /* is there a distinct solver at coarsest level? */
  if(sol!=NULL) {
    solver = sol;
  } else {
    solver = L[0] -> solver;
  }

  boil::oout << "Number of cycling levels: " << nlevels << boil::endl;
}
