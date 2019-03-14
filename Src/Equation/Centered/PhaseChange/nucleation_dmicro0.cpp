#include "nucleation.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
real Nucleation::dmicro0(const int i, const int j, const int k ) {
/***************************************************************************//**
*  \brief calculate initial thickness of micro layer
*  crude code: assume k-plane
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"pc.dmicro0: "<<boil::cart.iam()<<"\n";
#endif

  real rl=boil::exa;

  /* genuine sites */
  //for (int ns=0; ns < size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    if (sites[ns].active()) {
      real r = sqrt( pow(clr->xc(i)-sites[ns].x(),2.0)
                   + pow(clr->yc(j)-sites[ns].y(),2.0) );
      rl = min(rl,r);
    }
  }

  /* dummy sites */
  //for (int nsd=0; nsd < dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    if (dsites[nsd].active()) {
      real r = sqrt( pow(clr->xc(i)-dsites[nsd].x(),2.0)
                   + pow(clr->yc(j)-dsites[nsd].y(),2.0) );
      rl = min(rl,r);
    }
  }

  //real coef = 4.46e-3;  // Utaka's coefficient for water

  real d_return;
  if (rl==boil::exa) {
    d_return = boil::exa;
  } else if (rl < rmax) {
    //d_return = max(slope * rl, dmicro_min);
    d_return = max(slope * pow(rl, exp_slope), dmicro_min);
  } else {
    d_return = boil::exa;
  }

  return d_return;
}
