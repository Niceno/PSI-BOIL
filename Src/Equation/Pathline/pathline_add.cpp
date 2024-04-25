#include "pathline.h"

/***************************************************************************//**
*  add pathline
*******************************************************************************/
void Pathline::add(const real x, const real y, const real z) {
  id_serial++;
  //const Position pos(x,y,z);
  //Particle p(pos);
  Particle p(x, y, z, nval());
  //add((p));
  p.id(id_serial);
  particles.push_back(p);
  np(particles.size());

  if (np()%10==1)
  boil::oout<<"Pathline:add:np= "<<npa<<" x= "<<particles[npa-1].x()<<" y= "
            <<particles[npa-1].y()<<" z= "<<particles[npa-1].z()<<"\n";
}

