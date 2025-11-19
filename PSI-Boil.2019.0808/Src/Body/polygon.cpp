#include "body.h"

/******************************************************************************/
Polygon::Polygon(const int numb,
                 const real xp[], 
                 const real yp[], 
                 const real zp[], 
                 const real np[]) : nn(numb) {

  assert(numb >= 3); 
  assert(numb <= 6); 
  //assert(numb <= 12); 
 
  /* initialize */
  for_m(m) {
    g[~m]   = 0.0; 
    max[~m] = -FLT_MAX;          
    min[~m] =  FLT_MAX;      
  }

  /* copy node coordinates,
     gather center of gravity and
     bounding box coordinates */
  for(int n=0; n<nnodes(); n++) {
    x[n] = xp[n];
    y[n] = yp[n];
    z[n] = zp[n];
    g[0] += x[n];
    g[1] += y[n];
    g[2] += z[n];
    max[0] = boil::maxr(x[n], max[0]);
    min[0] = boil::minr(x[n], min[0]);
    max[1] = boil::maxr(y[n], max[1]);
    min[1] = boil::minr(y[n], min[1]);
    max[2] = boil::maxr(z[n], max[2]);
    min[2] = boil::minr(z[n], min[2]);
  }

  /* compute the center of gravity */
  for_m(m) g[~m] /= (real)nnodes();

  /* get normal */
  if(np) {
    for_m(m) nor[~m] = np[~m];

    this->arrange();

  /* compute normal */
  } else {
    real Ax = area_x();
    real Ay = area_y();
    real Az = area_z();
  
    real At = sqrt(Ax*Ax + Ay*Ay + Az*Az);

    //assert(At > 0.0);

    nor[0] = Ax/At;
    nor[1] = Ay/At;
    nor[2] = Az/At;
  }

  /* define equation for plane A x + B y + C z + D = 0 */
  real x10 = x[1]-x[0];  real y10 = y[1]-y[0];  real z10 = z[1]-z[0];
  real x20 = x[2]-x[0];  real y20 = y[2]-y[0];  real z20 = z[2]-z[0];
 
  A =   y10*z20 - y20*z10; 
  B = - x10*z20 + x20*z10; 
  C =   x10*y20 - x20*y10; 

  D = -x[0] * A - y[0] * B - z[0] * C;
}
