#include "heaviside.h"

/******************************************************************************/
void Heaviside::calculate(const bool evalflag) {
/***************************************************************************//**
*  \brief Calculate the Heaviside function from the color function.
*******************************************************************************/

  if(evalflag)
    evaluate_nodes();

  if(phi)
    calculate_heaviside(false);
  if(adens)
    calculate_adens(false);

  return;
}

void Heaviside::calculate_heaviside(const bool evalflag) {

  if(evalflag)
    evaluate_nodes();

  for_vijk((*phi),i,j,k) 
    (*phi)[i][j][k] = vf(i,j,k);
  
  (*phi).bnd_update();
  (*phi).exchange_all();

  return;
}

void Heaviside::calculate_adens(const bool evalflag) {

  if(evalflag)
    evaluate_nodes();

  for_vijk((*adens),i,j,k)
    (*adens)[i][j][k] = area(i,j,k)/adens->dV(i,j,k);
  (*adens).bnd_update();
  (*adens).exchange();
 
  return;
}
