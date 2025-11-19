updates
location.cpp
  introduce boil::BW
rack.h
  introduce 2D rack
Raviori/comp.x
  add Comp::x(), y(), z()
vof.cpp
  debug: nx->nz, mx->mz

main.cpp  (main-rack.cpp)

About Rack
 input range is 1-NX,1-NY,1-NZ (cell)
   output of scalar is cell
   output of vector is at face. uvw[Comp::u()][1][10][10]=0 if i-min is wall
 if the output is u-vector, i=NX+1 is theoretically possible,
   but it doesn't work due to implementation.

