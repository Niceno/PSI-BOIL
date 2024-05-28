#include "particle.h"
/******************************************************************************/
Particle::Particle(const real x, const real y, const real z, const int nval,
                   const real di, const real de) {
  xpos=x;
  ypos=y;
  zpos=z;
  dia = di;
  den = de;
  for (int ival = 0; ival < nval; ++ival) {
    sca.push_back(0.0);
  }
}
