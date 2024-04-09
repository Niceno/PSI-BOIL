#include "phasechangevof.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::finalize() {

#ifdef DEBUG
  boil::plot->plot(clr, tpr, phi, "clr-tpr-m",  time->current_step());
#endif

  /* calculate mdot */
  mdot();
#ifdef DEBUG
  boil::plot->plot(clr, tpr, phi, "clr-tpr-mdot",  time->current_step());
#endif

  /* calculate source terms */
  sources_clrs();
  sources_fext();
  sources_sum();
  boil::oout<<"phasechangevof_update: time= "<<time->current_time()
            <<" smdot_pos= "<<smdot_pos
            <<" smdot_neg= "<<smdot_neg<<"\n";
  smdot_pos_macro=smdot_pos;
  smdot_neg_macro=smdot_neg;

#ifdef DEBUG
  boil::plot->plot(clr, tprs, clrs, "clr-tprs-clrs",  time->current_step());
  std::cout<<"pc.update:end "<<boil::cart.iam()<<"\n";
  exit(0);
#endif

}

