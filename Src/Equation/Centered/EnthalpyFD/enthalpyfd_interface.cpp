#include "enthalpyfd.h"

/***************************************************************************//**
 *  Checks if the given cell is near an interface
******************************************************************************/

/*------------+
|  using vof  |
+------------*/
bool EnthalpyFD::Interface(const int dir, const Comp m,
                            const int i, const int j, const int k) {
  if(fs) {
    int of(0);
    if(dir>0)
      of = 1;
 
    if(m==Comp::i())
      return boil::realistic((*fs)[m][i+of][j][k]);  
    else if(m==Comp::j())
      return boil::realistic((*fs)[m][i][j+of][k]);  
    else
      return boil::realistic((*fs)[m][i][j][k+of]);  
  } else {
    int of(-1);
    if(dir>0)
      of = 1;
 
    real clrc = (*clr)[i][j][k];
    if(m==Comp::i())
      return ((clrc-clrsurf)*((*clr)[i+of][j][k]-clrsurf)<0.0);
    else if(m==Comp::j())
      return ((clrc-clrsurf)*((*clr)[i][j+of][k]-clrsurf)<0.0);
    else
      return ((clrc-clrsurf)*((*clr)[i][j][k+of]-clrsurf)<0.0);
  }

  return false;
}

bool EnthalpyFD::Interface_old(const int dir, const Comp m,
                                const int i, const int j, const int k) {
  if(fs) {
    int of(0);
    if(dir>0)
      of = 1;
 
    if(m==Comp::i())
      return boil::realistic(fsold[m][i+of][j][k]);  
    else if(m==Comp::j())
      return boil::realistic(fsold[m][i][j+of][k]);  
    else
      return boil::realistic(fsold[m][i][j][k+of]);  
  } else {
    int of(-1);
    if(dir>0)
      of = 1;
 
    real clrc = clrold[i][j][k];
    if(m==Comp::i())
      return ((clrc-clrsurf)*(clrold[i+of][j][k]-clrsurf)<0.0);
    else if(m==Comp::j())
      return ((clrc-clrsurf)*(clrold[i][j+of][k]-clrsurf)<0.0);
    else
      return ((clrc-clrsurf)*(clrold[i][j][k+of]-clrsurf)<0.0);
  }

  return false;
}
