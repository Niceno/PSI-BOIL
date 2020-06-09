#include "nucleation.h"

#if 0
/******************************************************************************/
real Nucleation::area_vapor(const Dir d,
                            const int i, const int j, const int k) const {
/******************************************************************************/
  Comp m;
  Sign sig;
  if       (d==Dir::imin()) {
    m = Comp::i();
    sig = Sign::neg();
  } else if(d==Dir::imax()) {
    m = Comp::i();
    sig = Sign::pos();
  } else if(d==Dir::jmin()) {
    m = Comp::j();
    sig = Sign::neg();
  } else if(d==Dir::jmax()) {
    m = Comp::j();
    sig = Sign::pos();
  } else if(d==Dir::kmin()) {
    m = Comp::k();
    sig = Sign::neg();
  } else if(d==Dir::kmax()) {
    m = Comp::k();
    sig = Sign::pos();
  }

  return area_vapor(sig,m,i,j,k);
}
#endif

/******************************************************************************/
real Nucleation::area_vapor(const Sign sig, const Comp & mcomp,
                            const int i, const int j, const int k) const {
/***************************************************************************//**
*  \brief approximate vapor area
*******************************************************************************/

  real frac = heavi->surface(Sign::neg(),Comp::k(),i,j,k);
  real ds = clr->dSz(Sign::neg(),i,j,k);

  if(sig==Sign::pos()) {
    return (1.0-frac) * ds;
  } else {
    return frac * ds;
  }
}
