#include "nucleation.h"

/***************************************************************************//**
*  give heat sink into wall
*******************************************************************************/
void Nucleation::heat_sink(std::vector<Site> & s, const int ns) {

  real vol_seed  = s[ns].vol_bubble();   // volume of seed
  real area_base = s[ns].area_base();  // area of bubble-base
  int kadj       = s[ns].kadj();          // k adjacent to wall
  //std::cout<<"vol_seed= "<<vol_seed<<" area_base= "<<area_base<<" "<<kadj<<"\n";

  real sum_q_sink = 0.0;  // positive sink, negative source

  //boil::oout<<"heat_sink:heat_sink()= "<<heat_sink()<<" "<<ns
  //          <<" contain "<<s[ns].contain_range()<<"\n";

  if(heat_sink()) {

    /* set qsrc */
    if (s[ns].contain_range()) {
      //boil::oout<<" i,j,k= "<<s[ns].is()<<" "<<s[ns].ie()
      //          <<" "<<s[ns].js()<<" "<<s[ns].je()<<" "<<kadj<<"\n";
      //boil::oout<<" vol_seed= "<<vol_seed<<" area_base= "<<area_base
      //          <<" "<<kadj<<"\n";
      for(int i=s[ns].is(); i<=s[ns].ie(); i++) {
        for(int j=s[ns].js(); j<=s[ns].je(); j++) {
          int k=kadj;
          real xcent=s[ns].x();
          real ycent=s[ns].y();
          real zcent=s[ns].z();
          real cseed = stratified_sphere(i,j,k,xcent,ycent,zcent);

          //if(matter_sig==Sign::neg()) {
          //  cseed = 1.-cseed;
          //}
          if (clr->domain()->ibody().off(i,j,k-1)) {
            if (area_base>0.0) {
              if (pre_heat_sink()) {
                /*  heat sink is given gradually in such a way that Tw
                    *  is always higher than Tsat */
                /* wall temperature */
		real Tw = (*tpr)[i][j][k-1];
		/* saturation temperature in adjacent cell */
		// real Tsat = cht->Tint(i,j,k);
		real Tsat = tsat;
                /* wall superheat */
                real wall_super = Tw - Tsat;
		/* heater power already calculated in main.cpp */
		/* THIS MEANS THAT nucleation.replant MUST BE CALLED IN MAIN.CPP
		 * AFTER DEFINING HEATER POWER
		 * BUT BEFORE SOLUTION OF ENERGY EQUATION */
		real q_heater = (*qsrc)[i][j][k-1]; // [W]
		/* power due to wall superheat */
                real p_wall_super = cps * wall_super * qsrc->dV(i,j,k-1) / time->dt(); //[W]
                /* heat sink in solid cell */
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
                                   / (area_base * time->dt() * qsrc->dzc(k-1)) 
                                   * (1.0-cseed) * qsrc->dV(i,j,k-1);
              }
            }
          }
        } /* loop per site cells */
      } /* loop per site cells */
    } /* contain_range */
  } /* is there qsrc */

  if (pre_heat_sink()) {
    /* energy of heat sink in this time step, in this decomposed domain */
    s[ns].set_sink_energy(sum_q_sink*time->dt());
    //std::cout<<"heat_sink: "<<s[ns].sink_energy()<<" [J], Energy for whole bubble "
    //         <<s[ns].vol_bubble()*rhov*latent<<"[J]\n";
  }

  return;
}
