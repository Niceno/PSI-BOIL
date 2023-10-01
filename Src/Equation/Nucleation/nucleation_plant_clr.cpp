#include "nucleation.h"

/***************************************************************************//**
*  plant single nucleation site
*******************************************************************************/
void Nucleation::plant_clr(const int ns) {

  if (sites[ns].contain_range()) {
    for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
      for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
        for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
          if(vf->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
            real xcent=sites[ns].x();
            real ycent=sites[ns].y();
            real zcent=sites[ns].z();

            real cseed = stratified_sphere(i,j,k, xcent,ycent,zcent);

            if(matter_sig==Sign::pos()) {
              (*vf)[i][j][k]=std::min((*vf)[i][j][k],cseed);
              (*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
            } else {
              (*vf)[i][j][k]=1.0-std::min(1.0-(*vf)[i][j][k],cseed);
              (*clr)[i][j][k]=1.0-std::min(1.0-(*clr)[i][j][k],cseed);
            }
          }
        }
      }
    }
  }
  //sites[ns].set_time_plant_clr(time->current_time()); // should not be called here!
  return;
}

