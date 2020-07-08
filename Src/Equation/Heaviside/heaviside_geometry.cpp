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

/***************************************************************************//**
*  \brief Cross product
*******************************************************************************/

/* calculate the cross product of the arguments */
Heaviside::XYZ Heaviside::CrossProduct(const XYZ & p1, const XYZ & p2) {
  return p1.CrossProduct(p2);
}

/***************************************************************************//**
*  \brief Shoelace formula (areapos is the centroid coords times area)
*******************************************************************************/

real Heaviside::Shoelace(const XY & v1, const XY & v2, const XY & v3,
                         XY & areapos) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v1.y - v1.x * v3.y);

  areapos.x = (v1.x * v2.y - v2.x * v1.y)*(v1.x+v2.x)
            + (v2.x * v3.y - v3.x * v2.y)*(v2.x+v3.x)
            + (v3.x * v1.y - v1.x * v3.y)*(v3.x+v1.x);

  areapos.y = (v1.x * v2.y - v2.x * v1.y)*(v1.y+v2.y)
            + (v2.x * v3.y - v3.x * v2.y)*(v2.y+v3.y)
            + (v3.x * v1.y - v1.x * v3.y)*(v3.y+v1.y);

  areapos.x /= 3.*area;
  areapos.y /= 3.*area;

  area = 0.5 * fabs(area);

  areapos.x *= area;
  areapos.y *= area;

  return(area);
}

real Heaviside::Shoelace(const XY & v1, const XY & v2, const XY & v3,
                         const XY & v4, XY & areapos) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v4.y - v4.x * v3.y)
       + (v4.x * v1.y - v1.x * v4.y);

  areapos.x = (v1.x * v2.y - v2.x * v1.y)*(v1.x+v2.x)
            + (v2.x * v3.y - v3.x * v2.y)*(v2.x+v3.x)
            + (v3.x * v4.y - v4.x * v3.y)*(v3.x+v4.x)
            + (v4.x * v1.y - v1.x * v4.y)*(v4.x+v1.x);

  areapos.y = (v1.x * v2.y - v2.x * v1.y)*(v1.y+v2.y)
            + (v2.x * v3.y - v3.x * v2.y)*(v2.y+v3.y)
            + (v3.x * v4.y - v4.x * v3.y)*(v3.y+v4.y)
            + (v4.x * v1.y - v1.x * v4.y)*(v4.y+v1.y);

  areapos.x /= 3.*area;
  areapos.y /= 3.*area;

  area = 0.5 * fabs(area);

  areapos.x *= area;
  areapos.y *= area;

  return(area);
}

real Heaviside::Shoelace(const XY & v1, const XY & v2, const XY & v3,
                         const XY & v4, const XY & v5, XY & areapos) {
  real area(0.0);

  area = (v1.x * v2.y - v2.x * v1.y)
       + (v2.x * v3.y - v3.x * v2.y)
       + (v3.x * v4.y - v4.x * v3.y)
       + (v4.x * v5.y - v5.x * v4.y)
       + (v5.x * v1.y - v1.x * v5.y);

  areapos.x = (v1.x * v2.y - v2.x * v1.y)*(v1.x+v2.x)
            + (v2.x * v3.y - v3.x * v2.y)*(v2.x+v3.x)
            + (v3.x * v4.y - v4.x * v3.y)*(v3.x+v4.x)
            + (v4.x * v5.y - v5.x * v4.y)*(v4.x+v5.x)
            + (v5.x * v1.y - v1.x * v5.y)*(v5.x+v1.x);

  areapos.y = (v1.x * v2.y - v2.x * v1.y)*(v1.y+v2.y)
            + (v2.x * v3.y - v3.x * v2.y)*(v2.y+v3.y)
            + (v3.x * v4.y - v4.x * v3.y)*(v3.y+v4.y)
            + (v4.x * v5.y - v5.x * v4.y)*(v4.y+v5.y)
            + (v5.x * v1.y - v1.x * v5.y)*(v5.y+v1.y);

  areapos.x /= 3.*area;
  areapos.y /= 3.*area;

  area = 0.5 * fabs(area);
  
  areapos.x *= area;
  areapos.y *= area;

  return(area);
}

/***************************************************************************//**
*  \brief Surface divergence
*******************************************************************************/

real Heaviside::triangle_surface_divergence(const TRIANGLE & t) {
  real sign(1.0);
  XYZ perpendicular = CrossProduct( t.p[1]-t.p[0] ,
                                    t.p[2]-t.p[0] );
  XYZ refvect[3];
    
  refvect[0] = t.v[0]-t.p[0]; 
  refvect[1] = t.v[1]-t.p[0]; 
  refvect[2] = t.v[2]-t.p[0]; 

  /* a given triangle has up to three reference points
     (grid vertices used for interpolating the position of the triangle
     which are below the isosurface); one can, however, be in fact
     located *above* or *on* the triangle plane as the surface can be
     concave. therefore, all three are checked and if less than two have
     the correct position, the normal vector sign is switched
     (i am not 100% sure this is a correct trigger, intution says yes) */

  int pos_count(0);
  if(perpendicular*refvect[0] < 0.0) /* correct */ 
    pos_count++;
  if(perpendicular*refvect[1] < 0.0) /* correct */ 
    pos_count++;
  if(perpendicular*refvect[2] < 0.0) /* correct */ 
    pos_count++;
  
  if(pos_count < 2)
    sign = -1.0;
    
  XYZ cross = CrossProduct(t.p[1],t.p[2]);
  real vol_area = 0.5 * sign * (cross*t.p[0]);

  return(vol_area);
}
