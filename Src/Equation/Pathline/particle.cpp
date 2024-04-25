#include "particle.h"

/******************************************************************************/
Particle::Particle(const real x, const real y, const real z, const int nval) {
  xpos=x;
  ypos=y;
  zpos=z;
  for (int ival = 0; ival < nval; ++ival) {
    sca.push_back(0.0);
  }
}

