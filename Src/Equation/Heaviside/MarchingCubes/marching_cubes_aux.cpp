#include "marching_cubes.h"
#include <iostream>

/***** 2D *****/

/* area calculation using the shoelace formula */

real MarchingCubes::SurfaceArea3(XY v1, XY v2, XY v3) {
  real area(0.0);
  
  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v1.y - v1.x * v3.y);

  area = 0.5 * fabs(area);

  return(area);
}

real MarchingCubes::SurfaceArea4(XY v1, XY v2, XY v3, XY v4) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v4.y - v4.x * v3.y)
       + (v4.x * v1.y - v1.x * v4.y);

  area = 0.5 * fabs(area);

  return(area);
}

real MarchingCubes::SurfaceArea5(XY v1, XY v2, XY v3, XY v4, XY v5) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v4.y - v4.x * v3.y)
       + (v4.x * v5.y - v5.x * v4.y)
       + (v5.x * v1.y - v1.x * v5.y);

  area = 0.5 * fabs(area);

  return(area);
}
