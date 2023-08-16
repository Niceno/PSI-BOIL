#include "pathline.h"

/***************************************************************************//**
*  add pathline
*******************************************************************************/
#if 0
void Pathline::add_local(const real x, const real y, const real z ) {
#else
void Pathline::add_local(const real x, const real y, const real z,
                         const real * di, const real * de) { 
#endif
#if 0
  if (di==NULL) {
    real dia = 0.0;
    real den = 0.0;
    particle_local.push_back({x, y, z, dia, den});
  } else {
    particle_local.push_back({x, y, z, * di, * de});
  }
#endif
  //std::cout<<"add_local "<<x<<" "<<y<<" "<<z<<"\n";
  //std::cout<<"add_local "<<di<<" "<<de<<"\n";   // 0 0
  //std::cout<<"add_local "<<&di<<" "<<&de<<"\n"; // 0x7fffffff6578 0x7fffffff6570
  //std::cout<<"add_local "<<*di<<" "<<*de<<"\n";
#if 0
  ParticleInit PI;
  PI.x=x;
  PI.y=y;
  PI.z=z;
  PI.diameter=di;
  PI.density=de;
  particle_local.push_back({PI});
#endif
  //particle_local.push_back({x, y, z, NULL, NULL});
  //*di = 0.0;
  //*de = 0.0;
  particle_local.push_back({x, y, z, di, de});
  //particle_local.push_back({x, y, z, di, de});

}

