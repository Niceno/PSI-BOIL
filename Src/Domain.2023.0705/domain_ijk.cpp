#include "domain.h"

/******************************************************************************/
int Domain::I(real x) const {

  for(int i_=0; i_<ni(); i_++) 
    if(x >= xn(i_) && x <= xn(i_+1)) return i_ + cxg().first() - 1;
  return -1;
}

/******************************************************************************/
int Domain::J(real y) const {

  for(int j_=0; j_<nj(); j_++) 
    if(y >= yn(j_) && y <= yn(j_+1)) return j_ + cyg().first() - 1;
  return -1;
}

/******************************************************************************/
int Domain::K(real z) const {

  for(int k_=0; k_<nk(); k_++) 
    if(z >= zn(k_) && z <= zn(k_+1)) return k_ + czg().first() - 1;
  return -1;
}
