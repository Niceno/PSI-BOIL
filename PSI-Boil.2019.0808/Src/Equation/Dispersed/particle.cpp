#include "particle.h"

/******************************************************************************/
Particle::Particle(const Position & pos, 
                   const Diameter & dia) {
  position = pos;
  velocity.set_to(0.0, 0.0, 0.0);
  diameter = dia;
  const real r = d() * 0.5;
  are = r * r * boil::pi;
  vol = 4.0 / 3.0 * r * r * r * boil::pi;

  boxdim.resize(DIM);
  boxuvw.resize(8);

  for(int c=0; c<8; c++) boxuvw[c].resize(DIM);
  voluvw.resize(DIM);
  volrho.resize(DIM);

  u_old.resize(DIM);
  uvwc_old(Comp::u(), 0.0);
  uvwc_old(Comp::v(), 0.0);
  uvwc_old(Comp::w(), 0.0);

  grad.resize(DIM);
  grad[0].resize(DIM);
  grad[1].resize(DIM);
  grad[2].resize(DIM);

}

/******************************************************************************/
Particle::Particle(const Position & pos, 
                   const Diameter & dia,  
                   const Position & vel) {
  position = pos;
  velocity = vel;                    
  diameter = dia;
  const real r = d() * 0.5;
  are = r * r * boil::pi;
  vol = 4.0 / 3.0 * r * r * r * boil::pi;

  u_old.resize(DIM);
  uvwc_old(Comp::u(), 0.0);
  uvwc_old(Comp::v(), 0.0);
  uvwc_old(Comp::w(), 0.0);

  grad.resize(DIM);
  grad[0].resize(DIM);
  grad[1].resize(DIM);
  grad[2].resize(DIM);

  boxdim.resize(DIM);
  boxuvw.resize(8);
  for(int c=0; c<8; c++) boxuvw[c].resize(DIM);
  voluvw.resize(DIM);
  volrho.resize(DIM);

}
