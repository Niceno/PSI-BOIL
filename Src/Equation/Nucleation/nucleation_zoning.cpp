#include "nucleation.h"

/******************************************************************************/
void Nucleation::zoning() {

  //boil::oout<<"nucleation:zoning\n";

  real xmin = vf->xn(1);
  real xmax = vf->xn(vf->ni());
  real lx = xmax-xmin;
  real ymin = vf->yn(1);
  real ymax = vf->yn(vf->nj());
  real ly = ymax-ymin;

  /* genuine sites */
  for (int ns=0; ns < size(); ns++){
    real xx = sites[ns].x();
    real yy = sites[ns].y();
    if(!limit_zoning ||
       ( (xmin-zoning_limit_multiplier*lx < xx) &&
         (xx < xmax+zoning_limit_multiplier*lx) &&
         (ymin-zoning_limit_multiplier*ly < yy) &&
         (yy < ymax+zoning_limit_multiplier*ly) 
       )
      ) {
      id_nearRegion.push_back(ns);
    }
  }
  //std::cout<<"zoning:id_nearRegion= "<<id_nearRegion.size()<<"\n";

  /* dummy sites */
  for (int nsd=0; nsd < dsize(); nsd++){
    real xx = dsites[nsd].x();
    real yy = dsites[nsd].y();
    if(!limit_zoning ||
       ( (xmin-zoning_limit_multiplier*lx < xx) &&
         (xx < xmax+zoning_limit_multiplier*lx) &&
         (ymin-zoning_limit_multiplier*ly < yy) &&
         (yy < ymax+zoning_limit_multiplier*ly)
       )
      ) {
      idd_nearRegion.push_back(nsd);
    }
  }
  //std::cout<<"zoning:idd_nearRegion= "<<idd_nearRegion.size()<<"\n";

  bzoning = true;
}
