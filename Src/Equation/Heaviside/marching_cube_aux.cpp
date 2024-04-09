#include "marching_cube.h"
#include <iostream>

/***** 2D *****/

/* linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value */
MarchingCube::XY MarchingCube::VertexInterp2D(real isolevel, XY p1, XY p2, real valp1, real valp2) {
  real mu;
  XY p;
  real eps;
  eps=1.0e-12;

  if (fabs(isolevel-valp1) < eps)
     return(p1);
  if (fabs(isolevel-valp2) < eps)
     return(p2);
  if (fabs(valp1-valp2) < eps)
     return(p1);

  mu = (isolevel - valp1) / (valp2 - valp1);
  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);

  return(p);
}

/* area calculation using the shoelace formula */

real MarchingCube::SurfaceArea3(XY v1, XY v2, XY v3) {
  real area(0.0);
  
  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v1.y - v1.x * v3.y);

  area = 0.5 * fabs(area);

  return(area);
}

real MarchingCube::SurfaceArea4(XY v1, XY v2, XY v3, XY v4) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v4.y - v4.x * v3.y)
       + (v4.x * v1.y - v1.x * v4.y);

  area = 0.5 * fabs(area);

  return(area);
}

real MarchingCube::SurfaceArea5(XY v1, XY v2, XY v3, XY v4, XY v5) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v4.y - v4.x * v3.y)
       + (v4.x * v5.y - v5.x * v4.y)
       + (v5.x * v1.y - v1.x * v5.y);

  area = 0.5 * fabs(area);

  return(area);
}

/***** 3D *****/

/* add two xyz vectors */
MarchingCube::XYZ MarchingCube::PlusXYZ(const XYZ p1, const XYZ p2) {
  XYZ pv;
  pv.x = p1.x + p2.x;
  pv.y = p1.y + p2.y;
  pv.z = p1.z + p2.z;

  return(pv);
}

/* returns the opposite xyz vector */
MarchingCube::XYZ MarchingCube::NegateXYZ(const XYZ p) {
  XYZ nv;
  nv.x = -p.x;
  nv.y = -p.y;
  nv.z = -p.z;

  return(nv);
}

/* calculate the dot product of the arguments */
real MarchingCube::DotProduct(const XYZ p1, const XYZ p2) {
  real dp(0.0);
  dp += p1.x * p2.x;
  dp += p1.y * p2.y;
  dp += p1.z * p2.z;

  return(dp);
}

/* calculate the cross product of the arguments */
MarchingCube::XYZ MarchingCube::CrossProduct(const XYZ p1, const XYZ p2) {
  XYZ cp;
  cp.x = p1.y * p2.z - p1.z * p2.y;
  cp.y = p1.z * p2.x - p1.x * p2.z;
  cp.z = p1.x * p2.y - p1.y * p2.x;

  return(cp);
}

/* linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value */
MarchingCube::XYZ MarchingCube::VertexInterp3D(real isolevel, XYZ p1, XYZ p2,
                                 real valp1, real valp2) {
   real mu;
   XYZ p;
   real eps;
   eps=1.0e-12;

   if (fabs(isolevel-valp1) < eps)
      return(p1);
   if (fabs(isolevel-valp2) < eps)
      return(p2);
   if (fabs(valp1-valp2) < eps)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p.x = p1.x + mu * (p2.x - p1.x);
   p.y = p1.y + mu * (p2.y - p1.y);
   p.z = p1.z + mu * (p2.z - p1.z);

   return(p);
}
