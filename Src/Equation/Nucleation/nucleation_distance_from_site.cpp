#include "nucleation.h"

/******************************************************************************/
real Nucleation::distance_from_site(const int i, const int j, const int k)
                                                                         const {
/***************************************************************************//**
*  \brief calculate distance from nucleation site
*  crude code: assume k-plane
*******************************************************************************/

  real rl = boil::unreal;

  /* genuine sites */
  //for (int ns=0; ns < size(); ns++){
#if 0
  if (i==26&&j==21&&k==15)
    std::cout<<"distance_from_site:26,21,15 "<<id_nearRegion.size()<<"\n";
#endif
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    if (sites[ns].active()) {
      real r = sqrt( pow(vf->xc(i)-sites[ns].x(),2.0)
                   + pow(vf->yc(j)-sites[ns].y(),2.0) );
      rl = std::min(rl,r);
    }
  }

  /* dummy sites */
  //for (int nsd=0; nsd < dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    if (dsites[nsd].active()) {
      real r = sqrt( pow(vf->xc(i)-dsites[nsd].x(),2.0)
                   + pow(vf->yc(j)-dsites[nsd].y(),2.0) );
      rl = std::min(rl,r);
    }
  }

  return rl;
}
