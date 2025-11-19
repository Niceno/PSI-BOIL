#include "nucleation.h"

/***************************************************************************//**
*  enforce to plant nucleation site regardless wall temperature etc.
*******************************************************************************/
void Nucleation::plant () {

  if (bzoning!=true) {
    zoning();
    //boil::oout<<"nucleation_plant:zoning\n";
  }

  real t_current = time->current_time();

  for(int ns=0; ns<size(); ns++){

    //boil::oout<<"nucl.plant: ID"<<ns<<" at: "<<sites[ns].x()<<", "
    //          <<sites[ns].y()<<", "<<sites[ns].z()<<"\n";

    /* time_Tact: time when Tw reached Tact (-> start heat sink) */
    sites[ns].set_time_Tact( t_current );
    /* set the site as active */
    sites[ns].set_active(true);

    if (!pre_heat_sink()) {
      /* if pre-heat-sink is not used, always set the flag false */
      sites[ns].set_qsink(false);

      /* if pre-heat-sink is not used, plant color function immidiately */
      plant_clr(sites,ns);
      sites[ns].set_plant_clr_current(true); /* flag for current time step */
      sites[ns].set_time_plant_clr( t_current ); /* store time */

      /* calculate heat sink */
      heat_sink(sites,ns);
      sites[ns].set_heat_sink_current(true); /* flag for current time step */

    } else {
      /* if pre-heat-sink is used, start giving heat-sink here.
       * color function will be planted in nucleation_replant. */

      /* calculate heat sink (clr is not modified in this function) */
      heat_sink(sites,ns);
      sites[ns].set_heat_sink_current(true); /* flag for current time step */
      /* set the flag active: if qsink()==true, heat_sink will be given
       * continuously */
      sites[ns].set_qsink(true);
    }
  } /* loop ns */

  /*---------------------------+
  |  dummy sites (periodic BC) |
  +---------------------------*/
  /* copy flags (plant_clr_current & heat_sink_current)
   * from father to children */

  for(int nsd=0; nsd<dsize(); nsd++){
  //for (int idd=0; idd<idd_nearRegion.size(); idd++){
    //int nsd=idd_nearRegion[idd];
    int fID=dsites[nsd].father();
    bool bclr = sites[fID].plant_clr_current();
    //boil::oout<<"plant:bclr0= "<<bclr<<" "<<fID<<" "<<idd<<"\n";
    dsites[nsd].set_plant_clr_current(bclr);
    //boil::oout<<"plant:bclr1= "<<bclr<<" "<<fID<<" "<<idd<<"\n";
    bool bhs = sites[fID].heat_sink_current();
    //boil::oout<<"plant:bclr,bhs0= "<<bclr<<" "<<bhs<<" "<<fID<<" "<<idd<<"\n";
    dsites[nsd].set_heat_sink_current(bhs);
    //boil::oout<<"plant:bclr,bhs1= "<<bclr<<" "<<bhs<<" "<<fID<<" "<<idd<<"\n";
  }

  /* plant color function */
  for(int nsd=0; nsd<dsize(); nsd++){
  //for (int idd=0; idd<idd_nearRegion.size(); idd++){
    //int nsd=idd_nearRegion[idd];
    if (dsites[nsd].plant_clr_current()) {
      plant_clr(dsites,nsd);
      //dsites[nsd].set_active(true);
    }
  }
  clr->exchange_all();

  /* heat sink */
  for(int nsd=0; nsd<dsize(); nsd++){
  //for (int idd=0; idd<idd_nearRegion.size(); idd++){
    //int nsd=idd_nearRegion[idd];
    boil::oout<<"dsites[nsd].heat_sink_current()= "<<dsites[nsd].heat_sink_current()<<"\n";
    if (dsites[nsd].heat_sink_current()) {
      heat_sink(dsites,nsd);
    }
  }
  qsrc->exchange();

  if (pre_heat_sink()) {
    /* sum sink_energy between procs.
     * sink_energy may be divided into sevral decomposed domains,
     * if the range covers several domains. */
    for (int ns=0; ns<size(); ns++){
      real tmp = sites[ns].sink_energy();
      boil::cart.sum_real(&tmp);
      sites[ns].set_sink_energy(tmp);
      //std::cout<<"plant:sink_energy= "<<tmp<<"\n";
    }
    /* sum energy in time */
    for (int ns=0; ns<size(); ns++){
      real tmp = sites[ns].sum_sink_energy();
      tmp += sites[ns].sink_energy();
      sites[ns].set_sum_sink_energy(tmp);
    }
  }

  /* set dmicro for new nucleation site */
  for(int ns=0; ns<size(); ns++){
  //for (int id=0; id<id_nearRegion.size(); id++){
    //int ns=id_nearRegion[id];
    //std::cout<<sites[ns].time_seed()<<" "<<t_current<<"\n";
    if( approx(sites[ns].time_plant_clr(),t_current,boil::pico) ){
      insert_dmicro(sites,ns);
    }
  }

  /* set dmicro for new dummy nucleation site */
  for(int nsd=0; nsd<dsize(); nsd++){
  //for (int idd=0; idd<idd_nearRegion.size(); idd++){
    //int nsd=idd_nearRegion[idd];
    int fID = dsites[nsd].father();
    if( approx(sites[fID].time_plant_clr(),t_current,boil::pico) ){
      insert_dmicro(dsites,nsd);
    }
  }

  deactivate_sites();

}
