#include "vof.h"

real VOF::extract_alpha(const int i, const int j, const int k) {

  /* step 0: is there interface? */
  real cc = phi[i][j][k];
  if(cc<boil::atto||cc>(1.0-boil::atto))
    return 0.0;

  /* step 1: calculate physical normal vector m0 */
  /* m0 points in the direction of decreasing vol. fraction */
  real dx = phi.dxc(i);
  real dy = phi.dyc(j);
  real dz = phi.dzc(k);

  real m0x = -nx[i][j][k]/dx;
  real m0y = -ny[i][j][k]/dy;
  real m0z = -nz[i][j][k]/dz;

  real m0norm = m0x*m0x+m0y*m0y+m0z*m0z;
  m0norm = sqrt(m0norm);

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

  real m2sum = m2x+m2y+m2z;

  m2x /= m2sum;
  m2y /= m2sum;
  m2z /= m2sum;
  
  /* step 4: calculate standardized alpha */
  real alp2 = calc_alpha(cc,m2x,m2y,m2z);

  /* step 5: destandardize the plane equation */
  real alp1 = alp2 * m2sum;

  /* step 6: remove the effect of mirroring */
  real alp0 = alp1;

  if(mirrorx) alp0 += m0x*dx; 
  if(mirrory) alp0 += m0y*dy; 
  if(mirrorz) alp0 += m0z*dz; 

  /* step 7: translate the coordinate system */
  real x0 = phi.xc(i) - 0.5*dx;
  real y0 = phi.yc(j) - 0.5*dy;
  real z0 = phi.zc(k) - 0.5*dz;

  alp0 += m0x*x0+m0y*y0+m0z*z0;

  /* the final equation has the form
   *  m0x*x + m0y*y + m0z*z = alp0
   *  (i.e. normal vector is pointing out of the liquid)
   */  
 
  return alp0;
}

