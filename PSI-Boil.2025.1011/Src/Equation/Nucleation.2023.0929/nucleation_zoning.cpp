#include "nucleation.h"

/******************************************************************************/
void Nucleation::zoning() {

  boil::oout<<"nucleation:zoning range_zoning="<<range_zoning<<"\n";

  real xmin = clr->xn(1);
  real xmax = clr->xn(clr->ni());
  real lx = xmax-xmin;
  real ymin = clr->yn(1);
  real ymax = clr->yn(clr->nj());
  real ly = ymax-ymin;

  const real range = range_zoning;  // 0.2 previously

  /* genuine sites */
  for (int ns=0; ns < size(); ns++){
    real xx = sites[ns].x();
    real yy = sites[ns].y();
    if (range>0) {
      if ((xmin-range*lx < xx) && (xx < xmax+range*lx)) {
      if ((ymin-range*ly < yy) && (yy < ymax+range*ly)) {
        id_nearRegion.push_back(ns);
      }
      }
    } else {
      id_nearRegion.push_back(ns);
    }
  }
  //std::cout<<"zoning: "<<boil::cart.iam()<<" "<<size()<<" "<<dsize()
  //         <<" "<<id_nearRegion.size()<<"\n";
  //exit(0);

  /* dummy sites */
  for (int nsd=0; nsd < dsize(); nsd++){
    real xx = dsites[nsd].x();
    real yy = dsites[nsd].y();
    if (range>0) {
      if ((xmin-range*lx < xx) && (xx < xmax+range*lx)) {
      if ((ymin-range*ly < yy) && (yy < ymax+range*ly)) {
        idd_nearRegion.push_back(nsd);
      }
      }
    } else {
      idd_nearRegion.push_back(nsd);
    }
  }
  //std::cout<<"zoning:idd_nearRegion= "<<idd_nearRegion.size()<<"\n";

  bzoning = true;
}
