#include "nucleation.h"

/***************************************************************************//**
*  plant single nucleation site
*******************************************************************************/
void Nucleation::plant_clr(std::vector<Site> & s, const int ns) {

  //boil::oout<<"plant_clr:begin "<<ns<<" "<<s[ns].contain_range()<<"\n";

  if (s[ns].contain_range()) {
    for(int i=s[ns].is(); i<=s[ns].ie(); i++) {
      for(int j=s[ns].js(); j<=s[ns].je(); j++) {
        for(int k=s[ns].ks(); k<=s[ns].ke(); k++) {
          if(clr->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
            real xcent=s[ns].x();
            real ycent=s[ns].y();
            real zcent=s[ns].z();

            real cseed = stratified_sphere(i,j,k, xcent,ycent,zcent);
	    if (fabs(cseed-clrsurf)<eps_clr) {
                cseed = clrsurf+copysign(1.0,cseed-clrsurf)*eps_clr;
            }

            //if(matter_sig==Sign::pos()) {
              (*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
            //} else {
            //  (*clr)[i][j][k]=1.0-std::min(1.0-(*clr)[i][j][k],cseed);
            //}
            if (cseed < threshold_c) (*tpr)[i][j][k]=tsat; // 2024.01.09
          }
        }
      }
    }
  }
  //sites[ns].set_time_plant_clr(time->current_time()); // should not be called here!

  return;
}

