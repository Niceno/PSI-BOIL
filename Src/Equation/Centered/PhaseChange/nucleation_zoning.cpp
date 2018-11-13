#include "nucleation.h"

/******************************************************************************/
void Nucleation::zoning() {

  //boil::oout<<"nucleation:zoning\n";

  real xmin = clr->xn(1);
  real xmax = clr->xn(clr->ni());
  real lx = xmax-xmin;
  real ymin = clr->yn(1);
  real ymax = clr->yn(clr->nj());
  real ly = ymax-ymin;

  /* genuine sites */
  for (int ns=0; ns < size(); ns++){
    real xx = sites[ns].x();
    real yy = sites[ns].y();
    if ((xmin-0.2*lx < xx) && (xx < xmax+0.2*lx)) {
    if ((ymin-0.2*ly < yy) && (yy < ymax+0.2*ly)) {
      id_nearRegion.push_back(ns);
    }
    }
  }
  //std::cout<<"zoning:id_nearRegion= "<<id_nearRegion.size()<<"\n";

  /* dummy sites */
  for (int nsd=0; nsd < dsize(); nsd++){
    real xx = dsites[nsd].x();
    real yy = dsites[nsd].y();
    if ((xmin-0.2*lx < xx) && (xx < xmax+0.2*lx)) {
    if ((ymin-0.2*ly < yy) && (yy < ymax+0.2*ly)) {
      idd_nearRegion.push_back(nsd);
    }
    }
  }
  //std::cout<<"zoning:idd_nearRegion= "<<idd_nearRegion.size()<<"\n";

  bzoning = true;
}

