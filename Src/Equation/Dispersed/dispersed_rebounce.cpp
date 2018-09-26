#include "dispersed.h"

/******************************************************************************/
void Dispersed::rebounce(int pa, int pb) {

  assert(pa > -1);
  assert(pa < size());
  assert(pb > -1);
  assert(pb < size());

  real dx = particles[pa].x() - particles[pb].x();
  real dy = particles[pa].y() - particles[pb].y();
  real dz = particles[pa].z() - particles[pb].z();
  real dist_ab  = sqrt(dx*dx + dy*dy + dz*dz);
  dx /= dist_ab; dy /= dist_ab; dz /= dist_ab;

  real ua_n = particles[pa].uvw(Comp::u()) * dx +   
              particles[pa].uvw(Comp::v()) * dy +     
              particles[pa].uvw(Comp::w()) * dz ;   
  real ub_n = particles[pb].uvw(Comp::u()) * dx +   
              particles[pb].uvw(Comp::v()) * dy +   
              particles[pb].uvw(Comp::w()) * dz ;   
 
  real mass_a = flu->rho(dispersed) * particles[pa].volume();
  real mass_b = flu->rho(dispersed) * particles[pb].volume();

  real delta_ua_n = 2.0 * (mass_a * ua_n + mass_b * ub_n) /
                          (mass_a + mass_b) - 2.0 * ua_n;    
  real delta_ub_n = 2.0 * (mass_a * ua_n + mass_b * ub_n) /
                          (mass_a + mass_b) - 2.0 * ub_n;    

  particles[pa].uvw(Comp::u()) += delta_ua_n * dx;    
  particles[pa].uvw(Comp::v()) += delta_ua_n * dy;    
  particles[pa].uvw(Comp::w()) += delta_ua_n * dz;    
  particles[pb].uvw(Comp::u()) += delta_ub_n * dx;    
  particles[pb].uvw(Comp::v()) += delta_ub_n * dy;    
  particles[pb].uvw(Comp::w()) += delta_ub_n * dz;   
}

/*-----------------------------------------------------------------------------+
 '$Id: dispersed_rebounce.cpp,v 1.1 2015/08/19 11:31:15 badreddine Exp $'/ 
+-----------------------------------------------------------------------------*/
