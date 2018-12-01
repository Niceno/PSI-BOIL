#include "enthalpytif.h"

/***************************************************************************//**
 *  Checks if the given cell is at an interface
******************************************************************************/

/*-------------------------+
| using indicator function |
+-------------------------*/
bool EnthalpyTIF::Interface(const int i, const int j, const int k,
                            const Scalar & heaviside) {
  if(Interface(heaviside[i][j][k])) {
    return true;
  } 
  return false;
}

bool EnthalpyTIF::Interface(const real heavi) {
  return heavi > boil::pico && heavi-1.0 < -boil::pico;
}

/*------------+
|  using vof  |
+------------*/
bool EnthalpyTIF::Interface(const int dir, const Comp m,
                            const int i, const int j, const int k) {
  int of(0);
  if(dir>0)
   of = 1;
 
  if(m==Comp::i())
    return (*fs)[m][i+of][j][k]<boil::zetta;  
  else if(m==Comp::j())
    return (*fs)[m][i][j+of][k]<boil::zetta;  
  else
    return (*fs)[m][i][j][k+of]<boil::zetta;  

  return false;
}
