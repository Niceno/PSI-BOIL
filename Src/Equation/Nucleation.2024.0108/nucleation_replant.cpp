#include "nucleation.h"

/***************************************************************************//**
*  plant nucleation site
*******************************************************************************/
void Nucleation::replant () {

  boil::timer.start("nucleation replant");

  if (bzoning!=true) {
    // Call only at initial time step or when restarted.
    // Zoning shortens CPU time for the calculation of
    // - distance from site to a cell (nucleation_distance_from_site.cpp),
    //   which is called from microlayer_d0.cpp in order to calculate 
    //   initial microlayer thickness.
    zoning();
  }

  real t_current = time->current_time();

  for(int ns=0; ns<size(); ns++){

    /* flag for current time step */
    sites[ns].set_plant_clr_current(false);  // default: plant color function
    sites[ns].set_heat_sink_current(false);  // default: plant color function

    if (pre_heat_sink()) {
      /* pre_heat_sink approach is used:
      *  Continue to give heat sink till sufficient energy is provided to bubble,
      *  then plant color function.
      *  Note that heat sink at each time step must be small enough to 
      *  satisfy Tw > Tsat. This is implemented in heat_sink */

      /* reset heat sink in this time step */
      sites[ns].set_sink_energy(0.0);
#if 0
      std::cout<<"replant:sites[ns].qsink()= "<<ns<<" "<<sites[ns].qsink()<<"\n";
#endif
      if (sites[ns].qsink()) {
        /* heat sink must be given, inside qsink period (see Powerpoint) */
        /* nucleation.pptx */
            
        /* calculate energy for the bubble */
        real energy_bubble = sites[ns].vol_bubble()*rhov*latent;

#if 0
        if (ns==0) {
          boil::oout<<"replant:ns="<<ns<<" energy_bubble= "<<energy_bubble
                    <<" sum_sink_energy= "<<sites[ns].sum_sink_energy()
                    <<" vol_bubble()= "<<sites[ns].vol_bubble()<<"\n";
        }
#endif
          
        if (sites[ns].sum_sink_energy()>=energy_bubble) {
          /* time-integrated heat-sink reaches bubble energy */

          boil::oout<<"replant:sum_sink_energy()>=energy_bubble ns="<<ns<<
                      " allow= "<<sites[ns].allow_replant()<<"\n";

          /* check if time step is small enough:
          *  = check answers from outside of class to allow replant */
          if ( sites[ns].allow_replant() ){
            /* plant color function */
            plant_clr(sites,ns);
            sites[ns].set_plant_clr_current(true);
            sites[ns].set_time_plant_clr( t_current );
            boil::oout<<"replant:plant_clr "<<t_current<<" ns "
                      <<ns<<" x "<<sites[ns].x()<<" y "<<sites[ns].y()<<"\n";

	    /* insert micro-layer */
	    insert_dmicro(sites,ns);

            /* reset sum_sink_energy after plant_clr */
            sites[ns].set_sum_sink_energy(0.0);

            // exit period for heat sink
            sites[ns].set_qsink(false);

            //boil::plot->plot(*cht->topo->clr,"clr",time->current_step());
          }

          /* request outside to allow replant */
          sites[ns].set_req_replant( true );
          //std::cout<<"replant: Pattern1 "<<boil::cart.iam()<<"\n";

        } else {
          /* heat sink < bubble energy */
          heat_sink(sites,ns);
          sites[ns].set_heat_sink_current(true); /* flag for current time step */
#if 0
          if(ns==4793){
            boil::oout<<"replant:sink_energy= "<<sites[ns].sink_energy()<<"\n";
          }
#endif
        }

      } else {
        /* outside of qsink period: two possibilities
	 * (1) bubble grows over the nucleation site
	 * (2) bubble does not exist and Tw < Tsat */

        //boil::oout<<"replant:active()= "<<ns<<" "<<sites[ns].active()<<"\n";
        if (sites[ns].active()) {
          /* bubble exists over the site */
        } else {
          /* bubble does not exist over the site
	   * check if need to start plant or not */

          /*---------------------------------------------------------------+
          |  start heat sink if (1) Tw > Tact,                             |
          |  and (2) seed point is liquid,                                 |
          |  and (3) bottom of previous bubble is higher than zplant()     |
          +---------------------------------------------------------------*/
          bool bseed=false;              /* bseed=true, if replant */
          real tpr_seed = tpr_site(ns);  /* Tw */
          real clr_seed = clr_site(ns);  /* color function at site */
          bool clr_seed_cond = !in_vapor(clr_seed); /* true: in liquid */
          bool bheight = height_bubble(ns);
          //boil::oout<<"tpr_seed= "<<tpr_seed<<" Tact "<<sites[ns].active_tpr()
	  //<<" bheight "<<bheight<<"\n";

          if( clr_seed_cond && tpr_seed > sites[ns].active_tpr() && bheight ) {
            /* three conditions (1-3) are satisfied */

            /* set flag active(): 
	     * true: give heat-sink or bubble exists over the site */
            sites[ns].set_active(true);
            boil::oout<<"replant:set_active:true<-false "<<t_current<<" ns "
                      <<ns<<" x "<<sites[ns].x()<<" y "<<sites[ns].y()<<"\n";

            /* time_Tact: time when Tw is supposed to reached Tact
	     *  (-> start giving heat sink) */
            sites[ns].set_time_Tact( t_current );

            /* if pre-heat-sink is used, set the flag qsink active 
	     * while heat-sink is given */
            sites[ns].set_qsink(true);
            /* calculate heat sink (clr is not modified in this function) */
            heat_sink(sites,ns);
            sites[ns].set_heat_sink_current(true); /* for current time step */
          }
        }
      }
#ifndef USE_VOF_NUCL
      /*--------------------------------------------------+ 
      |  continue replant during seed_period for CIPCSL2  |
      +--------------------------------------------------*/
      if( sites[ns].time_plant_clr() < t_current &&
          t_current < (sites[ns].time_plant_clr() + seed_period)) {
        plant_clr(sites,ns);
        sites[ns].set_plant_clr_current(true);  /* flag for current time step */

	/* set_time_plant_clr should not be call here */
	/* micro-layer must not be change here: don't call insert_dmicro */

        //boil::oout<<"replant:plant_clr "<<sites[ns].time_plant_clr()<<" "
        //          <<t_current<<" "<<seed_period<<"\n";
      }
#endif

    } else {
      /* pre_heat_sink approach is not used:
      *  plant color function and give heat sink in the same time step
      *  This may cause Tw < Tsat, if bubble energy is large
      *  e.g. high vapor density due to high pressure */

      if( sites[ns].time_plant_clr() < t_current &&
          t_current < (sites[ns].time_plant_clr() + seed_period)) {
        /*---------------------------------------+ 
        |  continue replant during seed_period   |
        +---------------------------------------*/
        plant_clr(sites,ns);
        sites[ns].set_plant_clr_current(true);
        if (size()==1) {
          boil::oout<<"replant:siteInfo "<<t_current<<" "<<tpr_site(ns)<<" "
                    <<clr_site(ns)<<" Pattern 1\n";
	}

	/* set_time_plant_clr should not be call here */
	/* micro-layer must not be change here: don't call insert_dmicro */

        /* request outside to allow replant */
        sites[ns].set_req_replant( true );

      } else {
        /* time outside of seed_period */
        /*---------------------------------------------------------------+
        |  If  (1) Tw > Tact,                                            |
        |  and (2) seed point is liquid,                                 |
        |  and (3) bottom of previous bubble is higher than zplant(),    |
        |  then plant color function and give heat sink                  |
        +---------------------------------------------------------------*/
        bool bseed=false;              /* bseed=true for replant */
        real tpr_seed = tpr_site(ns);  /* Tw */
        real clr_seed = clr_site(ns);  /* clr over the site */
        bool clr_seed_cond = !in_vapor(clr_seed); /* true: in liquid */
        bool bheight = height_bubble(ns); /* bubble bottom height */
        //boil::oout<<"tpr_seed= "<<tpr_seed<<" Tact "<<sites[ns].active_tpr()
	//<<" bheight "<<bheight<<"\n";
        if (size()==1) {
          boil::oout<<"replant:siteInfo "<<t_current<<" "<<tpr_seed<<" "
                    <<clr_seed<<" Pattern 2\n";
	}

        if( tpr_seed > sites[ns].active_tpr() && clr_seed_cond && bheight) {
          /* conditions (1), (2) and (3) are satisfied */

          if ( sites[ns].allow_replant() ){
            /* allow to replant (dt is small enough to plant) */
            //boil::oout<<"replant:Pattern2\n";

            /* set flag active() */
            sites[ns].set_active(true);
            boil::oout<<"replant:set_active:true<-false "<<t_current<<" ns "
                      <<ns<<" x "<<sites[ns].x()<<" y "<<sites[ns].y()<<"\n";

            /* plant color function */
	    plant_clr(sites,ns);
            sites[ns].set_plant_clr_current(true);
            sites[ns].set_time_plant_clr(t_current);

	    /* insert micro-layer */
	    insert_dmicro(sites,ns);

	    /* give heat sink */
            heat_sink(sites,ns);
            sites[ns].set_heat_sink_current(true); /* flag current time step */

            sites[ns].set_req_replant( true );
          } else {
            /* not allowed to replant (dt is not small enough to plant) */
            //boil::oout<<"replant:Pattern3, ns="<<ns<<"\n";
            sites[ns].set_req_replant( true );
          }
        } else {
          sites[ns].set_req_replant( false );
          //boil::oout<<"replant:Pattern4\n";
        }
      }
    }
  }

  /* for parallelization */
  clr->exchange_all();
  qsrc->exchange();

#if 0
  if (pre_heat_sink()) {
    /* sum sink_energy between procs.
     * sink_energy may be divided into sevral decomposed domains,
     * if the range covers several domains. */
    for (int ns=0; ns<size(); ns++){
      real tmp = sites[ns].sink_energy();
      boil::cart.sum_real(&tmp);
      sites[ns].set_sink_energy(tmp);
      //boil::oout<<"replant:sink_energy= "<<tmp<<"\n";
    }
    /* sum energy in time */
    for (int ns=0; ns<size(); ns++){
      real tmp = sites[ns].sum_sink_energy();
      tmp += sites[ns].sink_energy();
      sites[ns].set_sum_sink_energy(tmp);
      //boil::oout<<"replant:sum_sink_energy= "<<tmp<<"\n";
    }
  }
#endif

  /*--------------------------+
  |  copy bseed to dummy      |
  +--------------------------*/
  /* copy flags (plant_clr_current & heat_sink_current)
   * from father to children */

  for(int nsd=0; nsd<dsize(); nsd++){
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
    if (dsites[nsd].plant_clr_current()) {
      plant_clr(dsites,nsd);
    }
  }
  clr->exchange_all();

  /* heat sink */
  for(int nsd=0; nsd<dsize(); nsd++){
    dsites[nsd].set_sink_energy(0.0);
    if (dsites[nsd].heat_sink_current()) {
      heat_sink(dsites,nsd);  //dsites[nsd].set_sink_energy() is computed here
      //if (nsd==2) {std::cout<<"proc= "<<boil::cart.iam()<<" "<<dsites[nsd].sink_energy()<<"\n";}
    }
  }
  qsrc->exchange();
#if 0
  exit(0);
#endif

#if 0
  int ns=4793;
  boil::oout<<"replant:(1):ns="<<ns
              <<" sum_sink_energy= "<<sites[ns].sum_sink_energy()<<"\n";
#endif
  if (pre_heat_sink()) {
    /* sum sink_energy between procs.
     * sink_energy may be divided into sevral decomposed domains,
     * if the range covers several domains. */
    for (int ns=0; ns<size(); ns++){
      real tmp = sites[ns].sink_energy();
      boil::cart.sum_real(&tmp);
      sites[ns].set_sink_energy(tmp);
    }
    /* add sink_energy in dummy to that in genuine */
    for (int nsd=0; nsd<dsize(); nsd++) {
      int fID = dsites[nsd].father();
      real tmpd = dsites[nsd].sink_energy();
      //boil::oout<<"nsd= "<<nsd<<" fID= "<<fID<<" tmpd= "<<tmpd<<"\n";
      //if (nsd==2) {std::cout<<"proc= "<<boil::cart.iam()<<" "<<tmpd<<"\n";}
      boil::cart.sum_real(&tmpd);
      //boil::oout<<"nsd= "<<nsd<<" fID= "<<fID<<" tmpd= "<<tmpd<<"\n";
      real tmp = sites[fID].sink_energy();
      sites[fID].set_sink_energy(tmp+tmpd);
    }
#if 0
      int ns=4793;
      boil::oout<<"replant:(3):ns="<<ns
                <<" sink_energy= "<<sites[ns].sink_energy()<<"\n";
#endif
    /* sum energy in time */
    for (int ns=0; ns<size(); ns++){
      real tmp = sites[ns].sum_sink_energy();
      tmp += sites[ns].sink_energy();
      sites[ns].set_sum_sink_energy(tmp);
    }
  }

#if 0
  ns=4793;
  boil::oout<<"replant:(4):ns="<<ns
              <<" sum_sink_energy= "<<sites[ns].sum_sink_energy()<<"\n";
#endif
  /* set dmicro for new dummy nucleation site */
  for(int nsd=0; nsd<dsize(); nsd++){
    int fID = dsites[nsd].father();
    if( approx(sites[fID].time_plant_clr(),t_current,boil::pico) ){
      insert_dmicro(dsites,nsd);
    } 
  }

  deactivate_sites();

#if 0
  //int ns=4793;
  for(int ns=0; ns<sites.size(); ns++) {
    boil::oout<<"replant:end:ns= "<<ns<<" sum_sink_energy= "<<sites[ns].sum_sink_energy()<<"\n";
  }
#endif

  boil::timer.stop("nucleation replant");
}

