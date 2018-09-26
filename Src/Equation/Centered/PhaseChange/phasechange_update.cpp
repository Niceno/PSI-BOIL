#include "phasechange.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void PhaseChange::update(const Scalar * diff_eddy) {

  phi=0.0;
  /*------------------------------+
  |  calculate distance function  |
  +------------------------------*/
  /* calculate distance function */
  distfunc(clr,24);

  /* calculate normal vectors from distance function */
  gradphic(dist);
#ifdef DEBUG
  boil::plot->plot(clr, tpr, dist, "clr-tpr-dist",  time->current_step());
  boil::plot->plot(nx, ny, nz, "nx-ny-nz",  time->current_step());
#endif

  /* calculate grad(tpr) */
  gradt(diff_eddy);

  /* calculate m */
  m(diff_eddy);
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
  //sources_tprs();
  sources_sum();
  boil::oout<<"phasechange_update: time= "<<time->current_time()
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

/*-----------------------------------------------------------------------------+
 '$Id: phasechange_update.cpp,v 1.5 2015/05/05 14:48:35 sato Exp $'/
+-----------------------------------------------------------------------------*/
