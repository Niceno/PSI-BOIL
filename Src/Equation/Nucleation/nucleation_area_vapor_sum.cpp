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
  for (int i=vf->si(); i<=vf->ei(); i++) {
    if (vf->xc(i)<xr.first()) continue;
    if (vf->xc(i)>xr.last() ) continue;
    for (int j=vf->sj(); j<=vf->ej(); j++) {
      if (vf->yc(j)<yr.first()) continue;
      if (vf->yc(j)>yr.last() ) continue;
      for (int k=vf->sk(); k<=vf->ek(); k++) {
        if (vf->zc(k)<zr.first()) continue;
        if (vf->zc(k)>zr.last() ) continue;

        area += area_vapor(Sign::neg(),Comp::k(),i,j,k);
      }
    }
  }

  boil::cart.sum_real(&area);
  return(area);
}
