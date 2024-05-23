#include "nucleation.h"
#include "header.h"

/***************************************************************************//**
*  plant nucleation site
*******************************************************************************/
void Nucleation::replant () {

  //if (bzoning!=true) {
  if (bzoning==false) {
    // call only at initial time step or when restarted
    zoning();
  }
  //std::cout<<"replant:size= "<<id_nearRegion.size()<<"\n";

  /*-----------------+
  |  check criteria  |
  +-----------------*/

  real t_current = time->current_time();
  real front_array[6];

  for(int ns=0; ns<size(); ns++){

    bool bseed=false;              // bseed=true, if replant.
    real tpr_seed = tpr_site(ns);  // seed temperature
    real clr_seed = clr_site(ns);  // color function at seed point
    real zft;
    //boil::oout<<"nucleation_replant:ns= "<<ns<<" tpr_seed= "<<tpr_seed<<" "<<sites[ns].active_tpr()<<"\n";
    
    bool bheight = false;          // bheight=true, if zplant is satisfied.
    if(sites[ns].zplant()<0.0) {   // crude code
      bheight = true;
    } else {
      cht->topo->front_minmax(Range<real>(sites[ns].x()-dxmin,sites[ns].x()+dxmin),
                              Range<real>(sites[ns].y()-dxmin,sites[ns].y()+dxmin),
                              Range<real>(-boil::unreal,boil::unreal),
                              front_array);
      zft = front_array[4]; /* not optimal solution!!! */

      if(zft>sites[ns].zplant()){
        bheight=true;
      }
    }

    /*---------------------------------------------------------------+
    |  replant if (1) seed temp. is higher than activation temp.     |
    |  and (2) did not replant in the last time step                 |
    |  and (3) seed point is liquid  (time_seed affect cutneck)      |
    |  and (4) bottom of previous bubble is higher than zplant()     |
    |  and (5) avoid replant immediately after cutneck.              |
    +---------------------------------------------------------------*/
    bool clr_seed_cond = !in_vapor(clr_seed);
    //boil::oout<<"clr_seed,clr_seed_cond= "<<clr_seed<<" "<<clr_seed_cond<<"\n";
    if( tpr_seed > sites[ns].active_tpr() && 
        sites[ns].seed_prev()==false &&
        clr_seed_cond &&
        bheight &&
        //t_current > (sites[ns].time_cutneck() + period_prevent_replant) ) {
        t_current > (sites[ns].time_seed() + period_prevent_replant) ) {

      /* check answers from outside of class to allow replant */
      if ( sites[ns].allow_replant() ){
        sites[ns].set_time_seed( t_current );
        bseed=true;
        //boil::oout<<"replant:ns,x,y= "<<t_current<<" "<<ns<<" "
        //          <<sites[ns].x()<<" "<<sites[ns].y()<<"\n";
      }

      /* request outside to allow replant */
      sites[ns].set_req_replant( true );
      //std::cout<<"replant: Pattern1 "<<boil::cart.iam()<<"\n";

    } else {
      sites[ns].set_req_replant( false );
      //std::cout<<"replant: Pattern2 "<<boil::cart.iam()<<"\n";
    }

#ifndef USE_VOF_NUCL
    /*---------------------------------------+ 
    |  continue replant during seed_period   |
    +---------------------------------------*/
    if( t_current < (sites[ns].time_seed() + seed_period)) {
      bseed=true;
      sites[ns].set_req_replant( true );
      //std::cout<<"replant: Pattern2 "<<boil::cart.iam()<<"\n";
    }
#endif

    sites[ns].set_seed(bseed);
    //sites[ns].set_seed_prev(bseed);

    if (size()==1) {
      boil::oout<<"replant:printAll "<<ns<<" "<<t_current<<" "
        <<bseed<<" "<<tpr_seed<<" "<<sites[ns].active_tpr()<<" "
        <<clr_seed<<" "<<sites[ns].time_seed()<<" "
        <<sites[ns].allow_replant()<<" "<<sites[ns].req_replant()<<" "
        <<sites[ns].seed_prev()<<"\n"; 
    }
    //boil::oout<<"nucleation_replant:ns= "<<ns<<" bseed= "<<bseed<<"\n";
  }

  /*----------------------+
  |  copy bseed to dummy  |
  +----------------------*/
  for(int nsd=0; nsd<dsize(); nsd++){
    int ns=dsites[nsd].father();
    bool bseed = sites[ns].seed();
    dsites[nsd].set_seed(bseed);
  }

  /*----------+
  |  replant  |
  +----------*/
  /* genuine sites */
  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    if (sites[ns].seed() ) {
      sites[ns].set_active(true);
      plant_site(ns,!sites[ns].seed_prev());
    }
    sites[ns].set_seed_prev(sites[ns].seed());
  }

  /* dummy sites */
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    if (dsites[nsd].seed() ) {
      dsites[nsd].set_active(true);
      plant_dummy_site(nsd);
    }
  }

  vf->exchange_all();

  st_active();

#if 0
  boil::plot->plot(*cht->topo->clr, 
                   "clr-before",  time->current_step());
  exit(0);
#endif

  /* used by microlayer: first, update at walls must be called */
  // need to update heavi
  // Note: clr doesn't change due to heavi->calculate()
  heavi->calculate();
  //update_at_walls(conc.color(),conc.topo->fs);

#if 0
  boil::plot->plot(*cht->topo->clr, 
                   "clr-after",  time->current_step());
  exit(0); 
#endif

  // if microlayer is defined, dmicro will be set here
  // heavi->surface will be used to calculate area_vapor in upkeep_after_seeding
  //std::cout<<"replant:size= "<<id_nearRegion.size()<<"\n";
  //exit(0);
  upkeep_after_seeding();
#if 0
  boil::plot->plot(*cht->topo->clr, 
                   "clr-after",  time->current_step());
  exit(0);
#endif

  return;
}
