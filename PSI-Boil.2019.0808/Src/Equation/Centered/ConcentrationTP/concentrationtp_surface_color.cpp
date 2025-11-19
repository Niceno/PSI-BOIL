#include "concentrationtp.h"

/***************************************************************************//**
*  \brief Calculates the value of color function at the given surface.
*******************************************************************************/
real ConcentrationTP::surface_color(const Sign sig, const Comp & m,
                                    const int i, const int j, const int k) {
  return heavi->surface(sig,m,i,j,k);
}

real ConcentrationTP::surface_color(const Comp & m,
                                    const int i, const int j, const int k) {
  return heavi->surface(Sign::neg(),m,i,j,k);
}

