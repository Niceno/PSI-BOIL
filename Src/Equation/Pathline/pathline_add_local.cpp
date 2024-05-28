#include "pathline.h"

/***************************************************************************//**
*  add pathline
*******************************************************************************/
void Pathline::add_local(const real x, const real y, const real z,
                         const real dia, const real den ) {

  particle_local.push_back({x, y, z, dia, den});

}

