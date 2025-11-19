#include "nucleation.h"
using namespace std;

/******************************************************************************/
void Nucleation::optdat() {

  if (boptdat) {

    /* open a file */
    ofstpr.open("nucl-tpr.out",ios_base::app);
    ofsclr.open("nucl-clr.out",ios_base::app);

    if (boil::cart.iam()==0) {
      ofstpr<<time->current_step()<<" "<<time->current_time()<<" ";
      ofsclr<<time->current_step()<<" "<<time->current_time()<<" ";
    }

    for(int ns=0; ns<size(); ns++){
      real tpr_seed = tpr_site(ns);
      real clr_seed = clr_site(ns);
      if (boil::cart.iam()==0) {
        ofstpr<<tpr_seed<<" ";
        ofsclr<<clr_seed<<" ";
      }
    }

    if (boil::cart.iam()==0) {
      ofstpr<<"\n";
      ofsclr<<"\n";
    }

    /* close files */
    ofstpr.close();
    ofsclr.close();
  }
}

