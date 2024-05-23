#include "nucleation.h"

/***************************************************************************//**
*  Change the z-coordinate of nucleation sites.
*  This is necessary when you change radius of seed
*   (maybe due to change of grid resolution)
*******************************************************************************/
void Nucleation::change_site_z (real znew) {

  for(int ns=0; ns<size(); ns++){
    sites[ns].set_z(znew);
  } /* loop per sites */

  /* dummy sites (for symmetry) */
  for(int nsd=0; nsd<dsize(); nsd++){
    int ns=dsites[nsd].father();
    dsites[nsd].set_z(znew);
  }

  return;
}
