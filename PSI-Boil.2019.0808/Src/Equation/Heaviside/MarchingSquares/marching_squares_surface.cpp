#include "marching_squares.h"

/* calculate line fraction of cell face */
real MarchingSquares::surface(const Sign sig, const Comp & mcomp,
                              const int i, const int j, const int k) {
  int of(0);
  if(sig>0)
    of = 1;

  int im,ip,jm,jp,km,kp;
  /* degenerate case */
  if(mcomp==perpendicular) {
    return 0.0;
  } else {
    if       (perpendicular==Comp::i()) {
      im = i;
      ip = i;
      if(mcomp==Comp::j()) {
        jm = jp = j+of;
        km = k;
        kp = k+1;
      } else {
        km = kp = k+of;
        jm = j;
        jp = j+1;
      }
    } else if(perpendicular==Comp::j()) {
      jm = j;
      jp = j;
      if(mcomp==Comp::i()) {
        im = ip = i+of;
        km = k;
        kp = k+1;
      } else {
        km = kp = k+of;
        im = i;
        ip = i+1;
      }
    } else {
      km = k;
      kp = k;
      if(mcomp==Comp::i()) {
        im = ip = i+of;
        jm = j;
        jp = j+1;
      } else {
        jm = jp = j+of;
        im = i;
        ip = i+1;
      }
    }
  }

  real nodem = nodalvals[im][jm][km];
  real nodep = nodalvals[ip][jp][kp];

  if((nodem-clrsurf)*(nodep-clrsurf)>=0.0) {
    return nodem>clrsurf;
  } else {
    real node1 = std::min(nodem,nodep);
    real node2 = std::max(nodem,nodep);
    return (node2-clrsurf)/(node2-node1); 
  }

}
