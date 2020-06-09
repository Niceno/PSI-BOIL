#include "nucleation.h"

/***************************************************************************//**
*  plant nucleation site
*******************************************************************************/
void Nucleation::plant () {

  if (bzoning!=true) {
    zoning();
  }

  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];

    //boil::oout<<"nucl.plant: ID"<<ns<<" at: "<<sites[ns].x()<<", "
    //          <<sites[ns].y()<<", "<<sites[ns].z()<<"\n";

    sites[ns].set_time_seed( time->current_time() );
    sites[ns].set_active(true);

    plant_site(ns);
  } /* loop per sites */

  if(qsrc)
    qsrc->exchange();

  /* dummy sites */
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    //boil::oout<<"nucl.plantdummy: ID"<<nsd<<" at: "<<dsites[nsd].x()<<", "
    //          <<dsites[nsd].y()<<", "<<dsites[nsd].z()<<"\n";
    dsites[nsd].set_active(true);

    plant_dummy_site(nsd);
  } /* loop per dummy sites */

  vf->exchange_all();

  st_active();

  /* used by microlayer */
  upkeep_after_seeding();

  return;
}
