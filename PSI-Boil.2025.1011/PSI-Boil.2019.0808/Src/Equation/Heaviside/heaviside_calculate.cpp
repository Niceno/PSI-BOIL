#include "heaviside.h"

/******************************************************************************/
void Heaviside::calculate(const bool evalflag) {
/***************************************************************************//**
*  \brief Calculate the Heaviside function from the color function.
*******************************************************************************/

  if(evalflag)
    evaluate_nodes();

  calculate_flag(false);

  if(phi)
    calculate_vf(false);
  if(adens)
    calculate_adens(false);

  return;
}

void Heaviside::calculate_vf(const bool evalflag) {

  if(evalflag)
    evaluate_nodes();

  for_vijk((*phi),i,j,k) { 
    (*phi)[i][j][k] = vf(i,j,k);
  }
  
  (*phi).bnd_update();
  (*phi).exchange_all();

  return;
}

void Heaviside::calculate_adens(const bool evalflag) {

  if(evalflag)
    evaluate_nodes();

  for_vijk((*adens),i,j,k) {
    (*adens)[i][j][k] = ad(i,j,k);
  }
  (*adens).bnd_update();
  (*adens).exchange();
 
  return;
}

void Heaviside::calculate_flag(const bool evalflag) {

  if(evalflag)
    evaluate_nodes();

  for_vijk(flag,i,j,k) {
    flag[i][j][k] = status(i,j,k);
  }
  flag.bnd_update();
  flag.exchange();

  return;
}
