#include "nucleation.h"

/******************************************************************************/
void Nucleation::st_active() {
/***************************************************************************//**
*  \brief set nucleation site active or not
*******************************************************************************/
  for (int ns=0; ns<size(); ns++){
    //std::cout<<"st_active():clr_site= "<<clr_site(ns)<<"\n";

    bool clr_active_cond = in_vapor(clr_site(ns));

    if (clr_active_cond) {
      //sites[ns].set_active(true);
    } else {
      sites[ns].set_active(false);
    }
  }

  for (int nsd=0; nsd<dsize(); nsd++){
    int fID = dsites[nsd].father();
    bool bf = sites[fID].active();
    dsites[nsd].set_active(bf);
  }

  return;
}
