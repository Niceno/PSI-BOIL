#include "nucleation.h"

/******************************************************************************/
void Nucleation::deactivate_sites() {
/***************************************************************************//**
*  \brief change sites[ns].active() to false, once the site is covered with
*   liquid phase.
*   Parallelization OK.
*******************************************************************************/

  if(pre_heat_sink()) {
    for (int ns=0; ns<size(); ns++){
      if (!sites[ns].qsink() ) {
        bool clr_active_cond = in_vapor(clr_site(ns));
        if (clr_active_cond) {
          //sites[ns].set_active(true);
        } else {
          sites[ns].set_active(false);
        }
      }
    }
  } else {
    for (int ns=0; ns<size(); ns++){
      //std::cout<<"st_active():clr_site= "<<clr_site(ns)<<"\n";
      bool clr_active_cond = in_vapor(clr_site(ns));
      if (clr_active_cond) {
        //sites[ns].set_active(true);
      } else {
        sites[ns].set_active(false);
      }
    }
  }

  for (int nsd=0; nsd<dsize(); nsd++){
    int fID = dsites[nsd].father();
    bool bf = sites[fID].active();
    dsites[nsd].set_active(bf);
  }

  return;

}
