#include "nucleation.h"

/***************************************************************************//**
*  give heat sink into wall
*******************************************************************************/
void Nucleation::heat_sink(const int ns, const bool seed_source) {

  real vol_seed  = sites[ns].vol_bubble();   // volume of seed
  real area_base = sites[ns].area_base();  // area of bubble-base
  int kadj       = sites[ns].kadj();          // k adjacent to wall
  //std::cout<<"vol_seed= "<<vol_seed<<" area_base= "<<area_base<<" "<<kadj<<"\n";

  real sum_q_sink = 0.0;  // positive sink, negative source
  
  if(qsrc&&seed_source) {

    /* set qsrc */
    if (sites[ns].contain_range()) {
      for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
        for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
          int k=kadj;
          real xcent=sites[ns].x();
          real ycent=sites[ns].y();
          real zcent=sites[ns].z();
          real cseed = stratified_sphere(i,j,k,xcent,ycent,zcent);

          if(matter_sig==Sign::neg()) {
            cseed = 1.-cseed;
          }
          if (vf->domain()->ibody().off(i,j,k-1)) {
            if (area_base>0.0) {
              if (pre_heat_sink()) {
                /*  heat sink is given gradually in such a way that Tw
                    *  is always higher than Tsat */
                // wall temperature
		real Tw = cht->tmp()[i][j][k-1];
		// saturation temperature in adjacent cell
		real Tsat = cht->Tint(i,j,k);
                // wall superheat
                real wall_super = Tw - Tsat;
		// heater power given in main.cpp
		real q_heater = (*qsrc)[i][j][k-1]; // [W]
		// power due to wall superheat
                real p_wall_super = cps * wall_super * vf->dV(i,j,k-1) / time->dt(); //[W]
                // heat sink in solid cell
		real q_sink = (1.0-cseed) * (q_heater+p_wall_super);
		sum_q_sink += q_sink;
		(*qsrc)[i][j][k-1] -= q_sink;
                //std::cout<<"wall_super= "<<wall_super<<" Tw "<<Tw
                //         <<" Tsat "<<Tsat<<"\n";
                //std::cout<<"cseed= "<<cseed<<" clr= "<<(*clr)[i][j][k]
                //         <<" vf= "<<(*vf)[i][j][k]<<"\n";
              } else {
                /* heat sink is given in one time step */
                (*qsrc)[i][j][k-1] -= rhov * vol_seed * latent
                                   / (area_base * time->dt() * vf->dzc(k-1)) 
                                   * (1.0-cseed) * vf->dV(i,j,k-1);
              }
            }
          }
        } /* loop per site cells */
      } /* loop per site cells */
    } /* contain_range */
  } /* is there qsrc */
  //exit(0);

  if (pre_heat_sink()) {
    /* energy of heat sink in this time step, in this decomposed domain */
    sites[ns].set_sink_energy(sum_q_sink*time->dt());
    //std::cout<<"heat_sink: "<<sites[ns].sink_energy()<<" [J], Energy for whole bubble "
    //         <<sites[ns].vol_bubble()*rhov*latent<<"[J]\n";
  }

  return;
}

