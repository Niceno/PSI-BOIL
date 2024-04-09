#include "nucleation.h"
using namespace std;

/******************************************************************************/
real Nucleation::area_vapor_sum( Range<real> xr
                               , Range<real> yr
                               , Range<real> zr ) {
/***************************************************************************//**
*  \brief calculate vapor area on wall
*           assumption: wall is flat
*******************************************************************************/

  //bool bib=true;  // true: with IB,  false: without IB
  //if ( clr->domain()->ibody().ncall()==0 )   bib=false;

  real c, ds, area;
  area = 0.0;
  for (int i=(*clr).si(); i<=(*clr).ei(); i++) {
    if ((*clr).xc(i)<xr.first()) continue;
    if ((*clr).xc(i)>xr.last() ) continue;
    for (int j=(*clr).sj(); j<=(*clr).ej(); j++) {
      if ((*clr).yc(j)<yr.first()) continue;
      if ((*clr).yc(j)>yr.last() ) continue;
      for (int k=(*clr).sk(); k<=(*clr).ek(); k++) {
        if ((*clr).zc(k)<zr.first()) continue;
        if ((*clr).zc(k)>zr.last() ) continue;

        c  = max(min((*clr)[i][j][k],1.0),0.0);
        ds = clr->dSz(i,j,k);
        area += (1.0-c) * ds;
	//std::cout<<i<<" "<<j<<" "<<k<<" "<<c<<" "<<ds<<" "
        //<<(*clr).zc(k)<<" "<<zr.first()<<" "<<zr.last()<<"\n";

      }
    }
  }
   boil::cart.sum_real(&area);
   return(area);
}
