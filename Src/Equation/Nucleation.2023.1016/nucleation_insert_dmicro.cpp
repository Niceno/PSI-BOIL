#include "nucleation.h"

/***************************************************************************//**
*  insert dmicro after plant of color function
*******************************************************************************/
void Nucleation::insert_dmicro(std::vector<Site> & s, const int ns) {

  if (s[ns].contain_range()) {
    for(int i=s[ns].is(); i<=s[ns].ie(); i++) {
      for(int j=s[ns].js(); j<=s[ns].je(); j++) {
        int k = s[ns].kadj(); // in fluid
        if (area_vapor(i,j,k,Dir::kmin()) >0.0 ) {
          //if(clr_new-clr_old<0.0){// decrease-of-liquid (=increase-of-vapor)
            // modify dmicro
            dmicro[i][j][k]=std::min(dmicro0(i,j,k),dmicro[i][j][k]);  //update on 2021.12.18
          //}
        }
      }
    }
  }
  //s[ns].set_time_plant_clr(time->current_time()); // should not be called here!
  return;
}

