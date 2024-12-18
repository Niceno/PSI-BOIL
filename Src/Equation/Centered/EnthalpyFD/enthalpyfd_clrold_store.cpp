#include "enthalpyfd.h"

/***************************************************************************//**
*  store clrold
*  clrold: color function in previous time step
*******************************************************************************/
void EnthalpyFD::clrold_store() {

  /* initial time step or restart */
  if(!store_clrold){
    boil::oout<<"EnthalpyFD::new_time_step()  initialize clrold"<<"\n";
    for_aijk(i,j,k){
      clrold[i][j][k] = (*clr)[i][j][k];
    }
    store_clrold = true;
  }

}
