#include "nucleation.h"
using namespace std;

/******************************************************************************/
void Nucleation::st_active() {
/***************************************************************************//**
*  \brief set nucleation site active or not
*******************************************************************************/
  for (int ns=0; ns<size(); ns++){
    //std::cout<<"st_active():clr_site= "<<clr_site(ns)<<"\n";
    if ( clr_site(ns) < 0.5) {
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

}
