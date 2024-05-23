#include "nucleation.h"

/***************************************************************************//**
*  Enforce to plant bubble regardless activation temperature nor time increment
*  If there are multiple nucleation sites, all of them will be activated
*******************************************************************************/
void Nucleation::plant () {

  if (bzoning!=true) {
    zoning();
  }

  for(int ns=0; ns<size(); ns++){

    //boil::oout<<"nucl.plant: ID"<<ns<<" at: "<<sites[ns].x()<<", "
    //          <<sites[ns].y()<<", "<<sites[ns].z()<<"\n";

    // time_Tact: time when Tw is supposed to reached Tact (-> start heat sink)
    sites[ns].set_time_Tact( time->current_time() );
    // mark the site as active
    sites[ns].set_active(true);

    if (!pre_heat_sink()) {
      // if pre-heat-sink is not used, set the flag false
      sites[ns].set_qsink(false);

      // if pre-heat-sink is not used, plant color function immidiately
      sites[ns].set_time_plant_clr( time->current_time() );
      plant_clr(ns);
      // calculate heat sink
      heat_sink(ns);
    } else {
      // if pre-heat-sink is used, set the flag active
      sites[ns].set_qsink(true);
      // calculate heat sink (clr is not modified in this function)
      heat_sink(ns);
    }
  } /* loop per sites */

  vf->exchange_all();
  clr->exchange_all();

  if(qsrc){
    qsrc->exchange();

    if (pre_heat_sink()) {
      /* sum sink_energy between procs.
       * sink_energy may be divided into sevral decomposed domains,
       * if the range covers several domains. */
      for (int ns=0; ns<size(); ns++){
        real tmp = sites[ns].sink_energy();
        boil::cart.sum_real(&tmp);
        sites[ns].set_sink_energy(tmp);
        std::cout<<"plant:sink_energy= "<<tmp<<"\n";
      }
      /* sum energy in time */
      for (int ns=0; ns<size(); ns++){
        real tmp = sites[ns].sum_sink_energy();
        tmp += sites[ns].sink_energy();
        sites[ns].set_sum_sink_energy(tmp);
      }
    }
  }

  // necessary for activate flag
  deactivate_sites();

  /* used by microlayer: first, update at walls must be called */
  //upkeep_after_seeding();  // this is called at the beginning of
                             // microlayer_update.cpp

  return;
}
