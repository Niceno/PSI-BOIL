#include "nucleation.h"
#include "header.h"
//#define DEBUG

/***************************************************************************//**
*  plant nucleation site
*******************************************************************************/
void Nucleation::replant () {

#ifdef DEBUG
  boil::oout<<"replant:begin\n";
#endif
  if (bzoning==false) {
    // Call only at initial time step or when restarted.
    // Zoning shortens CPU time for the calculation of
    // - distance from site to a cell (nucleation_distance_from_site.cpp),
    //   which is called from microlayer_d0.cpp in order to calculate 
    //   initial microlayer thickness.
    zoning();
  }
  //std::cout<<"replant:size= "<<id_nearRegion.size()<<"\n";

  real t_current = time->current_time();
  real front_array[6];

  for(int ns=0; ns<size(); ns++){
  //for (int id=0; id<id_nearRegion.size(); id++){  // don't use this
    //int ns=id_nearRegion[id];

    if (pre_heat_sink()) {
      // pre_heat_sink approach is used:
      // Continue to give heat sink till sufficient energy is provided to bubble,
      // then plant color function.
      // Note that heat sink at each time step must be small enough to 
      // satisfy Tw > Tsat. This is implemented in heat_sink

      // reset heat sink in this time step
      sites[ns].set_sink_energy(0.0);

      //boil::oout<<"replant:sites[ns].qsink()= "<<ns<<" "<<sites[ns].qsink()<<"\n";

      if (sites[ns].qsink()) {
        // heat sink must be given

        // update qsink
        real energy_bubble = sites[ns].vol_bubble()*rhov*latent;

        if (sites[ns].sum_sink_energy()>=energy_bubble) {
	  // heat sink reached bubble energy

          boil::oout<<"replant:sum_sink_energy()>=energy_bubble ns="<<ns<<"\n";
	  boil::oout<<"replant:allow_replant= "<<sites[ns].allow_replant()<<"\n";
          // check if time step is small enough
          /* check answers from outside of class to allow replant */
            if ( sites[ns].allow_replant() ){
              // plant color function
              plant_clr(ns);
              sites[ns].set_time_plant_clr( t_current );
              boil::oout<<"replant:plant_clr\n";

              // reset sum_sink_energy after plant_clr
              sites[ns].set_sum_sink_energy(0.0);

              // exit period for heat sink
              sites[ns].set_qsink(false);

              //boil::plot->plot(*cht->topo->clr,"clr",time->current_step());
            }

            /* request outside to allow replant */
            sites[ns].set_req_replant( true );
            //std::cout<<"replant: Pattern1 "<<boil::cart.iam()<<"\n";
        
        } else {
          // heat sink < bubble energy
          heat_sink(ns);
	}

      } else {
        // outside of qsink period

        //boil::oout<<"replant:active()= "<<ns<<" "<<sites[ns].active()<<"\n";
        if (sites[ns].active()) {
          // bubble exists over the site
        } else {
          // bubble does not exist over the site
	  // check if need to start plant or not

          bool bseed=false;              // bseed=true, if replant.
          real tpr_seed = tpr_site(ns);  // seed temperature
          real clr_seed = clr_site(ns);  // color function at seed point
          //boil::oout<<"tpr_seed= "<<tpr_seed<<" Tact "<<sites[ns].active_tpr()<<"\n";

          /*---------------------------------------------------------------+
          |  start heat sink if (1) Tw > Tact,                             |
          |  and (2) seed point is liquid,                                 |
          |  and (3) bottom of previous bubble is higher than zplant()     |
          +---------------------------------------------------------------*/
          bool clr_seed_cond = !in_vapor(clr_seed);

          bool bheight = false;          // bheight=true, if zplant is satisfied.
          if(sites[ns].zplant()<0.0) {   // crude code
            bheight = true;
          } else {
            cht->topo->front_minmax(Range<real>(sites[ns].x()-dxmin,sites[ns].x()+dxmin),
                              Range<real>(sites[ns].y()-dxmin,sites[ns].y()+dxmin),
                              Range<real>(-boil::unreal,boil::unreal),
                              front_array);
            real zft = front_array[4]; /* not optimal solution!!! */
            if(zft>sites[ns].zplant()){
              bheight=true;
            }
          }

          if( clr_seed_cond && tpr_seed > sites[ns].active_tpr() && bheight ) {

            // set flag active()
            sites[ns].set_active(true);
            boil::oout<<"replant:set_active(true) from active()==false ns= "<<ns<<"\n";

            // time_Tact: time when Tw is supposed to reached Tact (-> start heat sink)
            sites[ns].set_time_Tact( t_current );
            //boil::aout<<"replant:set_time_Tact= "<<sites[ns].time_Tact()<<"\n";

            // if pre-heat-sink is used, set the flag active
            sites[ns].set_qsink(true);
            // calculate heat sink (clr is not modified in this function)
            heat_sink(ns);
          }
        }
      }
#ifndef USE_VOF_NUCL
      /*---------------------------------------+ 
      |  continue replant during seed_period   |
      +---------------------------------------*/
      if( sites[ns].time_plant_clr() < t_current &&
          t_current < (sites[ns].time_plant_clr() + seed_period)) {
        plant_clr(ns);
        //boil::oout<<"replant:plant_clr "<<sites[ns].time_plant_clr()<<" "
        //          <<t_current<<" "<<seed_period<<"\n";
      }
#endif

    } else {
      // pre_heat_sink approach is not used:
      // plant color function and give heat sink in the same time step
      // This may cause Tw < Tsat, if bubble energy is large,
      // i.e. high vapor density due to high pressure.

      if( sites[ns].time_plant_clr() < t_current &&
          t_current < (sites[ns].time_plant_clr() + seed_period)) {
        /*---------------------------------------+ 
        |  continue replant during seed_period   |
        +---------------------------------------*/
        plant_clr(ns);
        boil::oout<<"replant:plant_clr:Pattern1 "<<sites[ns].time_plant_clr()<<" "
                  <<t_current<<" "<<seed_period<<"\n";

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
        bool bseed=false;              // bseed=true, if replant.
        real tpr_seed = tpr_site(ns);  // seed temperature
        real clr_seed = clr_site(ns);  // color function at seed point
        //boil::oout<<"tpr_seed= "<<tpr_seed<<" Tact "<<sites[ns].active_tpr()<<"\n";

        bool clr_seed_cond = !in_vapor(clr_seed);

        bool bheight = false;          // bheight=true, if zplant is satisfied.
        if(sites[ns].zplant()<0.0) {   // crude code
          bheight = true;
        } else {
          cht->topo->front_minmax(Range<real>(sites[ns].x()-dxmin,sites[ns].x()+dxmin),
                            Range<real>(sites[ns].y()-dxmin,sites[ns].y()+dxmin),
                            Range<real>(-boil::unreal,boil::unreal),
                            front_array);
          real zft = front_array[4]; /* not optimal solution!!! */
          if(zft>sites[ns].zplant()){
            bheight=true;
          }
        }

        if( tpr_seed > sites[ns].active_tpr() && clr_seed_cond && bheight) {
          /* conditions (1), (2) and (3) are satisfied */

          if ( sites[ns].allow_replant() ){
            /* allow to replant (dt is small enough to plant) */
	    boil::oout<<"replant:Pattern2\n";

            // set flag active()
            sites[ns].set_active(true);
            boil::oout<<"replant:set_active(true) from active()==false ns= "<<ns<<"\n";

            heat_sink(ns);  // give heat sink in one time step
            plant_clr(ns);  // plant color function

            sites[ns].set_time_plant_clr(t_current);
            sites[ns].set_req_replant( true );
	  } else {
            /* not allowed to replant (dt is not small enough to plant) */
	    boil::oout<<"replant:Pattern3\n";
            sites[ns].set_req_replant( true );
	  }
        } else {
          sites[ns].set_req_replant( false );
	  boil::oout<<"replant:Pattern4\n";
	}
      }
    }
  }

  /* for parallelization */
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
  }

  // necessary for activate flag
  deactivate_sites();
  //boil::oout<<"replant:time= "<<t_current<<" active "<<sites[0].active()
  //          <<" qsink "<<sites[0].qsink()
  //          <<" sum_sink_energy "<<sites[0].sum_sink_energy()
  //          <<" tpr "<<tpr_site(0)<<" clr "<<clr_site(0)<<"\n";

  upkeep_after_seeding();

#if 0
  if (size()==1) {
      int ns = 0;
      boil::oout<<"replant:printAll "<<ns<<" "<<t_current<<" "
        <<bseed<<" "<<tpr_seed<<" "<<sites[ns].active_tpr()<<" "
        <<clr_seed<<" "<<sites[ns].time_seed()<<" "
        <<sites[ns].allow_replant()<<" "<<sites[ns].req_replant()<<" "
        <<sites[ns].seed_prev()<<"\n";
  }
#endif
#ifdef DEBUG
  boil::oout<<"replant:end\n";
#endif
  return;
}
