#include "phasechangevof.h"

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
bool PhaseChangeVOF::Interface(const int i, const int j, const int k,
                               const Comp m, const int dir) {
  int of(1);
  if(dir<0) of = 0;

  if(m == Comp::i())
    return (fs)[m][i+of][j][k] < boil::zetta;
  else if(m == Comp::j()) 
    return (fs)[m][i][j+of][k] < boil::zetta; 
  else
    return (fs)[m][i][j][k+of] < boil::zetta; 

  return false;
}
