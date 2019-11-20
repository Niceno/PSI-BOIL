#include "enthalpyfd.h"

/***************************************************************************//**
 *  Checks if the given cell is near an interface
******************************************************************************/
bool EnthalpyFD::Interface(const int dir, const Comp m,
                            const int i, const int j, const int k) {
  int of(1);
  if(dir<0) of = 0;

  if(m==Comp::i())
    return boil::realistic(fs[m][i+of][j][k]);
  else if(m==Comp::j())
    return boil::realistic(fs[m][i][j+of][k]); 
  else
    return boil::realistic(fs[m][i][j][k+of]); 

  return false;
}

bool EnthalpyFD::Interface_old(const int dir, const Comp m,
                                const int i, const int j, const int k) {
  int of(1);
  if(dir<0) of = 0;
 
  if(m==Comp::i())
    return boil::realistic(fsold[m][i+of][j][k]);  
  else if(m==Comp::j())
    return boil::realistic(fsold[m][i][j+of][k]);  
  else
    return boil::realistic(fsold[m][i][j][k+of]);  

  return false;
}
