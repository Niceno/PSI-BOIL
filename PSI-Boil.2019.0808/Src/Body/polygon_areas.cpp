#include "body.h"

/******************************************************************************/
real Polygon::area_x() const {

  real area = 0.0;

  for(int i=0; i<nnodes()-1; i++) 
    area += (z[i]+z[i+1]) * (y[i]-y[i+1]);

  area += (z[nnodes()-1]+z[0]) * (y[nnodes()-1]-y[0]);

  return (area*0.5);
}

/******************************************************************************/
real Polygon::area_y() const {

  real area = 0.0;

  for(int i=0; i<nnodes()-1; i++) 
    area += (x[i]+x[i+1]) * (z[i]-z[i+1]);

  area += (x[nnodes()-1]+x[0]) * (z[nnodes()-1]-z[0]);

  return (area*0.5);
}

/******************************************************************************/
real Polygon::area_z() const {

  real area = 0.0;

  for(int i=0; i<nnodes()-1; i++) 
    area += (y[i]+y[i+1]) * (x[i]-x[i+1]);

  area += (y[nnodes()-1]+y[0]) * (x[nnodes()-1]-x[0]);

  return (area*0.5);
}
