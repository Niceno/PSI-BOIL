#include "additive.h"

/***************************************************************************//**
*  Constructs the AC solver. As a first step, coarsend the Centered variable
*  for which it is meant to be used. As a second step, stores pointers to
*  variable's levels.  
*******************************************************************************/
AC::AC(Centered * cen) {

  /* set variables which steer the v-cycle */
  max_cyc = 20;
  min_cyc =  0;
  targ_res_val = boil::nano;
  targ_res_rat = 0.01;
  stop_if_div  = true;

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

  boil::oout << "Number of cycling levels: " << nlevels << boil::endl;
}

/*-----------------------------------------------------------------------------+
 '$Id: additive.cpp,v 1.16 2011/09/22 13:19:19 niceno Exp $'/
+-----------------------------------------------------------------------------*/
