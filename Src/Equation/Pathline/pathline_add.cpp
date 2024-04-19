#include "pathline.h"

/***************************************************************************//**
*  add pathline
*******************************************************************************/
void Pathline::add(const real x, const real y, const real z) {
  //const Position pos(x,y,z);
  //Particle p(pos);
  Particle p(x, y, z);
  //add((p));
  particles.push_back(p);
  npa = particles.size(); /* take the label of the last particle */
  boil::oout<<"Pathline:add:np= "<<npa<<" x= "<<particles[npa-1].x()<<" y= "
            <<particles[npa-1].y()<<" z= "<<particles[npa-1].z()<<"\n";
}

