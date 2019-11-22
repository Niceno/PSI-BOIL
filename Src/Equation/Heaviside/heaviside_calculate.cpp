#include "heaviside.h"

/******************************************************************************/
void Heaviside::calculate() {
/***************************************************************************//**
*  \brief Calculate the Heaviside function from the color function.
*******************************************************************************/

  if(phi)
    calculate_heaviside();
  if(adens)
    calculate_adens();

  return;
}

void Heaviside::calculate_heaviside() {

  for_vijk((*phi),i,j,k) 
    (*phi)[i][j][k] = volume(i,j,k);
  
  (*phi).bnd_update();
  (*phi).exchange_all();

  return;
}

void Heaviside::calculate_adens() {

  for_vijk((*adens),i,j,k)
    (*adens)[i][j][k] = area(i,j,k)/adens->dV(i,j,k);
  (*adens).bnd_update();
  (*adens).exchange();
 
  return;
}
