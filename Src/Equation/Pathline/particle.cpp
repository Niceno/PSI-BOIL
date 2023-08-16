#include "particle.h"
/******************************************************************************/
Particle::Particle(const real x, const real y, const real z, const int nval,
                   const real * di, const real * de) {
  //std::cout<<"particle "<<x<<" "<<y<<" "<<z<<"\n";
  //std::cout<<x<<" "<<y<<" "<<z<<" "<<*di<<" "<<*de<<"\n";
  xpos=x;
  ypos=y;
  zpos=z;
  dia = 0.0;
  if(di!=NULL) dia = * di;
  //std::cout<<"particle "<<x<<" "<<y<<" "<<z<<"\n";
  den = 0.0;
  if(de!=NULL) den = *de;
  //den = * de;
  for (int ival = 0; ival < nval; ++ival) {
    sca.push_back(0.0);
  }
}
