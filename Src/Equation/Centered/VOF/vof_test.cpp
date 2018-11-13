#include "vof.h"

real VOF::test(const real DX, const real DY, const real DZ,
               const real MX, const real MY, const real MZ,
               const real X0, const real Y0, const real Z0,
               const real CC) {

  real mx = -MX;
  real my = -MY;
  real mz = -MZ;
  real dx = DX;
  real dy = DY;
  real dz = DZ;
  real x0 = X0;
  real y0 = Y0;
  real z0 = Z0;
  real cc = CC;

  real nx = mx*dx;
  real ny = my*dy;
  real nz = mz*dz;

  real m0x = mx;
  real m0y = my;
  real m0z = mz;

  real m0norm = m0x*m0x+m0y*m0y+m0z*m0z;
  m0norm = sqrt(m0norm);

#if 0
  boil::oout <<"VOF-test, m0norm= "<<m0norm<<boil::endl;
#endif 

  m0x /= m0norm;
  m0y /= m0norm;
  m0z /= m0norm;

  /* step 2: transform m0 to positive values by mirroring */
  bool mirrorx(false), mirrory(false), mirrorz(false);
  
  real m1x(m0x);
  real m1y(m0y);
  real m1z(m0z);

  if(m0x<0.0) { m1x = -m1x; mirrorx = true; }
  if(m0y<0.0) { m1y = -m1y; mirrory = true; }
  if(m0z<0.0) { m1z = -m1z; mirrorz = true; }

  /* step 3: transform m1 to the standardized space */
  real m2x(m1x*dx);
  real m2y(m1y*dy);
  real m2z(m1z*dz);


#if 1
  real m2sum = m2x+m2y+m2z;

  #if 0
    boil::oout <<"VOF-test, m2sum= "<<m2sum<<boil::endl;
  #endif 

  m2x /= m2sum;
  m2y /= m2sum;
  m2z /= m2sum;
#endif  

#if 0
  boil::oout <<"VOF-test, m2= "<<m2x<<" "<<m2y<<" "<<m2z<<boil::endl;
#endif 

  /* step 4: calculate standardized alpha */
  real alp2 = calc_alpha(cc,m2x,m2y,m2z);

#if 0
  boil::oout <<"VOF-test, alp2= "<<alp2<<boil::endl;
#endif 

  /* step 5: destandardize the plane equation */
  real alp1 = alp2 * m2sum;

#if 0
  boil::oout <<"VOF-test, alp1= "<<alp2<<boil::endl;
#endif 

  /* step 6: remove the effect of mirroring */
  real alp0 = alp1;

#if 1
  if(mirrorx) alp0 += m0x*dx; 
  if(mirrory) alp0 += m0y*dy; 
  if(mirrorz) alp0 += m0z*dz; 
#endif

#if 0
  boil::oout <<"VOF-test, alp0= "<<alp2<<boil::endl;
#endif 
 
  /* step 7: translate the coordinate system */
  alp0 += mx*x0+my*y0+mz*z0; 

  return alp0;
}

