#include "phasechange4.h"
//#define DEBUG

/******************************************************************************/
void PhaseChange4::finalize() {

#ifdef DEBUG
  boil::plot->plot(*(topo->clr), tpr, phi, "clr-tpr-m",  time->current_step());
#endif

  /* calculate mdot */
  mdot();
#ifdef DEBUG
  boil::plot->plot(*(topo->clr), tpr, phi, "clr-tpr-mdot",  time->current_step());
#endif

  /* calculate source terms */
  sources();
  boil::oout<<"phasechange4_update: time= "<<time->current_time()
            <<" smdot_pos= "<<smdot_pos
            <<" smdot_neg= "<<smdot_neg<<"\n";

#ifdef DEBUG
  boil::plot->plot(*(topo->clr), tprs, vfs, "clr-tprs-vfs",  time->current_step());
  std::cout<<"pc.update:end "<<boil::cart.iam()<<"\n";
  exit(0);
#endif

}

