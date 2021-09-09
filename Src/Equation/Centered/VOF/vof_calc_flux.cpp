#include "vof.h"

real VOF::calc_flux(const real g, real c, 
                    const real nx, const real ny, const real nz,
                    const real nalpha) {

  /* no advection */
  if(g==0.0) {
    return 0.0;
  }

  /* n points to the liquid */
  real vn1 = -nx;
  real vn2 = -ny;
  real vn3 = -nz;

  /* cell cut */
  real absg = fabs(g);

  /* L1-normalized normal vector */
  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3)+boil::pico;
  real qa = 1.0/(vm1+vm2+vm3);
  vm1 *= qa;
  vm2 *= qa;
  vm3 *= qa;
  
  /* normalized alpha value in the cell
   * even if nalpha is not realistic, the checks at the beginning
   * of calc_v should handle it */
  real alpha = qa*nalpha;
      
  /* this number serves two purposes:
     - alpha offset due to origin shift
     - 1-ra is the new L1 norm of the vm vector */
  real ra = vm1 * (1.0 - absg);
  qa = 1.0/(1.0-ra);

  /* if directions are aligned, alpha offset has to be introduced
     if not, the mirroring of the coordinate system nullifies the offset  */
  if (g*vn1 > 0) alpha = alpha -ra;

  /* normal vector in direction of transport is rescaled due to
     augmentation of cell dimensions */
  vm1 = vm1 * absg;

  /* calculate f: reduced flux (for unit volume)
     - the qa factor renormalizes the line equation
     - the g factor then scales the transported volume back to the org cube 
     - it also imposes the correct sign */
  real f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;

  return f;
}

