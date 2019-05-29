#include "domain.h"

/******************************************************************************/
bool Domain::bnd_symmetry(const Dir d) const {
  if(d==Dir::imin()) {
    return imins;
  } else if(d==Dir::imax()) {
    return imaxs;
  } else if(d==Dir::jmin()) {
    return jmins;
  } else if(d==Dir::jmax()) {
    return jmaxs;
  } else if(d==Dir::kmin()) {
    return kmins;
  } else if(d==Dir::kmax()) {
    return kmaxs;
  }
}
