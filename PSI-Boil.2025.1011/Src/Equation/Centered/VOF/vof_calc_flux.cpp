#include "vof.h"

real VOF::calc_flux(const real g, real c, const real nx, const real ny, const real nz) {

  if(g==0.0) {
    return 0.0;
  }

  real vn1 = -nx;
  real vn2 = -ny;
  real vn3 = -nz;

  real absg = fabs(g);
  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3)+boil::pico;
  real qa = 1.0/(vm1+vm2+vm3);
  vm1 *= qa;
  vm2 *= qa;
  vm3 *= qa;
  
  real alpha = calc_alpha(c, vm1, vm2, vm3);

  //std::cout<<alpha<<" ";
      
  real ra = vm1 * (1.0 - absg);
  qa = 1.0/(1.0-ra);
  if (g*vn1 > 0) alpha = alpha -ra;
  vm1 = vm1 * absg;

  /* calculate f: reduced flux (for unit volume) */
  real f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;

  return f;
}

