#include "phasechangevof.h"

/******************************************************************************/
real PhaseChangeVOF::distance_center(const Sign sig, const Comp & m,
                                     const int i, const int j, const int k) {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell center  
*******************************************************************************/
  if       (m==Comp::i()) {
    if(sig<0) {
      return phi.dxw(i);
    } else {
      return phi.dxe(i);
    }
  } else if(m==Comp::j()) {
    if(sig<0) {
      return phi.dys(j);
    } else {
      return phi.dyn(j);
    }
  } else {
    if(sig<0) {
      return phi.dzb(k);
    } else {
      return phi.dzt(k);
    }
  }

  return 0.0;
}
