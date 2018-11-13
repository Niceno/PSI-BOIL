#include "heaviside.h"

/******************************************************************************/
void Heaviside::calculate() {
/***************************************************************************//**
*  \brief Calculate the Heaviside function from the color function.
*******************************************************************************/

  for_vijk(phi,i,j,k) 
    phi[i][j][k] = mc.volume(i,j,k);
  
  phi.bnd_update();
  phi.exchange_all();

  return;
}
