#include "domain.h"

/******************************************************************************/
bool Domain::bnd_symmetry(const Dir d) const {

  int axis(0), dir(0);

  if       (d==Dir::imin()) {
    axis = 0;
    dir  = 0;
  } else if(d==Dir::imax()) {
    axis = 0;
    dir  = 1;
  } else if(d==Dir::jmin()) {
    axis = 1;
    dir  = 0;
  } else if(d==Dir::jmax()) {
    axis = 1;
    dir  = 1;
  } else if(d==Dir::kmin()) {
    axis = 2;
    dir  = 0;
  } else if(d==Dir::kmax()) {
    axis = 2;
    dir  = 1;
  }

  return ctf[axis][dir];
}
