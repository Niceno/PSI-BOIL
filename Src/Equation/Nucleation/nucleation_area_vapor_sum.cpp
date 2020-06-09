#include "nucleation.h"

/******************************************************************************/
real Nucleation::area_vapor_sum( Range<real> xr
                               , Range<real> yr
                               , Range<real> zr ) {
/***************************************************************************//**
*  \brief approximate vapor area on wall
*           assumption: wall is flat, in the negative z-direction
*******************************************************************************/

  real area(0.0);
  for (int i=(*clr).si(); i<=(*clr).ei(); i++) {
    if ((*clr).xc(i)<xr.first()) continue;
    if ((*clr).xc(i)>xr.last() ) continue;
    for (int j=(*clr).sj(); j<=(*clr).ej(); j++) {
      if ((*clr).yc(j)<yr.first()) continue;
      if ((*clr).yc(j)>yr.last() ) continue;
      for (int k=(*clr).sk(); k<=(*clr).ek(); k++) {
        if ((*clr).zc(k)<zr.first()) continue;
        if ((*clr).zc(k)>zr.last() ) continue;

        area += area_vapor(Sign::neg(),Comp::k(),i,j,k);
      }
    }
  }

  boil::cart.sum_real(&area);
  return(area);
}
