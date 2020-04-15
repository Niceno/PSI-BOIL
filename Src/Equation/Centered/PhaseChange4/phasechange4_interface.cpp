#include "phasechange4.h"

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
bool PhaseChange4::interface(const Sign dir, const Comp m,
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

bool PhaseChange4::interface(const int i, const int j, const int k) {
  if(adens[i][j][k]>boil::pico) {
    return true;
  }

  return false;
}

