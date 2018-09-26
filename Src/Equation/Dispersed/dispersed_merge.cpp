#include "dispersed.h"

/******************************************************************************/
void Dispersed::merge(int p) {

  assert(p > -1);
  assert(p < size());

  /* browse in particle "p" box */
  for_pijk(p, i, j, k) {

    const real c = interface_fraction(i, j, k, p, continuous);
    if(dispersed_fraction(c) > 0.0) (*col)[i][j][k] = c;
    
  }
   erase(p);
   (*col).exchange_all();
}

/******************************************************************************/
void Dispersed::merge(int pa, int pb) {

  assert(pa > -1);
  assert(pa < size());
  assert(pb > -1);
  assert(pb < size());

  /* compute volume after merging */
  const real merged_volume = particles[pa].volume() 
                           + particles[pb].volume(); 

  /* compute weighted coordinates and velocity after merging */
  real merged_xyz[DIM] = {0.0, 0.0, 0.0};
  real merged_uvw[DIM] = {0.0, 0.0, 0.0};

  for_m(m) {
    merged_xyz[~m] += ( particles[pa].volume() * particles[pa].xyz(m) +
                        particles[pb].volume() * particles[pb].xyz(m) ) 
                      / merged_volume; 
    merged_uvw[~m] += ( particles[pa].volume() * particles[pa].uvw(m) +
                        particles[pb].volume() * particles[pb].uvw(m) ) 
                      / merged_volume; 
  }

  /* compute diameter after merging, delete old and add a new particle */
  const real merged_diameter = 2.0 * pow((.75*merged_volume/boil::pi),(1./3.));
  OPR(merged_diameter); 
 
  /* must erase higher index first.  
     otherwise, size decreases, and 
     higher index is out of range */
  erase( boil::maxi(pa,pb) );
  erase( boil::mini(pa,pb) );

  /* add new particle */
  particles.push_back( 
    Particle( 
      Position(merged_xyz[0], merged_xyz[1], merged_xyz[2]),
      Diameter(merged_diameter),
      Position(merged_uvw[0], merged_uvw[1], merged_uvw[2]) 
    ) 
  );
   
}

/*-----------------------------------------------------------------------------+
 '$Id: dispersed_merge.cpp,v 1.4 2015/08/19 11:39:12 badreddine Exp $'/ 
+-----------------------------------------------------------------------------*/
