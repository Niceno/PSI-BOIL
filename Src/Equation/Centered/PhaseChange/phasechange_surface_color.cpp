#include "phasechange.h"

/***************************************************************************//**
*  \brief Calculates the value of color function at the given surface.
*******************************************************************************/
real PhaseChange::surface_color(const Comp & m,
                                const int i, const int j, const int k) {
  real ds;
  if (m == Comp::u())
    ds = (*bndclr).dSx(m,i,j,k);
  else if (m == Comp::v())
    ds = (*bndclr).dSy(m,i,j,k);
  else
    ds = (*bndclr).dSz(m,i,j,k);

  return (*bndclr)[m][i][j][k] / ds;
}

