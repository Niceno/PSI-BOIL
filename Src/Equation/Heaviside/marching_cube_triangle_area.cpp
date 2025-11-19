#include "marching_cube.h"

/* calculate area contribution from a triangular surface */
real MarchingCube::triangle_area(const TRIANGLE t) {
  real x1 = t.p[1].x-t.p[0].x;
  real y1 = t.p[1].y-t.p[0].y;
  real z1 = t.p[1].z-t.p[0].z;
  real x2 = t.p[2].x-t.p[0].x;
  real y2 = t.p[2].y-t.p[0].y;
  real z2 = t.p[2].z-t.p[0].z;
  real area=  (y1*z2-z1*y2)*(y1*z2-z1*y2)
             +(z1*x2-x1*z2)*(z1*x2-x1*z2)
             +(x1*y2-y1*x2)*(x1*y2-y1*x2);
  area = 0.5 * sqrt(area);
  return(area);
}

/* calculate volumetric contribution from a triangular surface */
real MarchingCube::triangle_vol_area(const TRIANGLE t) {
  real sign(1.0);
  XYZ perpendicular = CrossProduct( PlusXYZ(t.p[1] , NegateXYZ(t.p[0])) ,
                                    PlusXYZ(t.p[2] , NegateXYZ(t.p[0])) );
  XYZ refvect[3];
    
  refvect[0] = PlusXYZ(t.v[0],NegateXYZ(t.p[0])); 
  refvect[1] = PlusXYZ(t.v[1],NegateXYZ(t.p[0])); 
  refvect[2] = PlusXYZ(t.v[2],NegateXYZ(t.p[0])); 

  /* a given triangle has up to three reference points
     (grid vertices used for interpolating the position of the triangle
     which are below the isosurface); one can, however, be in fact
     located *above* or *on* the triangle plane as the surface can be
     concave. therefore, all three are checked and if less than two have
     the correct position, the normal vector sign is switched
     (i am not 100% sure this is a correct trigger, intution says yes) */

  int pos_count(0);
  if (DotProduct(perpendicular,refvect[0]) < 0.0) /* correct */ 
       pos_count++;
  if (DotProduct(perpendicular,refvect[1]) < 0.0) /* correct */ 
       pos_count++;
  if (DotProduct(perpendicular,refvect[2]) < 0.0) /* correct */ 
       pos_count++;
  
  if (pos_count < 2) sign = -1.0;
    
  XYZ cross = CrossProduct(t.p[1],t.p[2]);
  real vol_area = 0.5 * sign * DotProduct(cross,t.p[0]);

  return(vol_area);
}
