#include "heaviside.h"

/***************************************************************************//**
*  \brief Vertex interpolation
*******************************************************************************/

/* linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value */
Heaviside::XY Heaviside::VertexInterp(const real & isolevel, 
                                      const XY & p1, const XY & p2,
                                      const real & valp1, const real & valp2) {

  real mu;
  XY p;

  if (fabs(isolevel-valp1) < boil::pico)
     return(p1);
  if (fabs(isolevel-valp2) < boil::pico)
     return(p2);
  if (fabs(valp1-valp2) < boil::pico)
     return(p1);

  mu = (isolevel - valp1) / (valp2 - valp1);
  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);

  return(p);
}

Heaviside::XYZ Heaviside::VertexInterp(const real & isolevel,
                                       const XYZ & p1, const XYZ & p2,
                                       const real & valp1, const real & valp2) {
   real mu;
   XYZ p;

   if (fabs(isolevel-valp1) < boil::pico)
      return(p1);
   if (fabs(isolevel-valp2) < boil::pico)
      return(p2);
   if (fabs(valp1-valp2) < boil::pico)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p.x = p1.x + mu * (p2.x - p1.x);
   p.y = p1.y + mu * (p2.y - p1.y);
   p.z = p1.z + mu * (p2.z - p1.z);

   return(p);
}

/* cross product */
/* calculate the cross product of the arguments */
Heaviside::XYZ Heaviside::CrossProduct(const XYZ & p1, const XYZ & p2) {
  return p1.CrossProduct(p2);
}
