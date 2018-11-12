#include "cipcsl2.h"
#include <iomanip>
using namespace std;

/******************************************************************************/
void CIPCSL2::redist(const bool local) {
/***************************************************************************//**
*  \brief Reset color function using re-distance technique and/or hyperbolic
*         tangent function.
*         input & output: clr
*         tempolary: sclr
*******************************************************************************/

#ifdef DEBUG
  boil::oout<<"cipcsl2_redist:start\n";
#endif

  /*--------------------------------+
  |  check necessity of redistance  |
  +--------------------------------*/
  /* calculate distance function */
  //distfunc(clr,24);  // distfunc is calculated in cipcsl2_advance()

  /* calculate normal vectors from distance function */
  gradphic(dist);

  /* calculate alp. alp will be used in sharpen(). */
  set_alp();
  //  boil::plot->plot((*u),clr,alp, "uvw-clr-alp", time->current_step());

  /* copy clr to sclr */
  for_aijk(i,j,k)
    sclr[i][j][k]=clr[i][j][k];

  /*-------------+
  |  sharpening  |
  +-------------*/
  //std::cout.setf(std::ios_base::scientific);
  //std::cout<< setprecision(16);
  //boil::oout<<"before "<<totalvol(clr)<<"\n";
  if(sharpen(sclr, 0.5, itsharpen, local)){
    /* update clr and scheme variables */
    update_node(sclr);
  }
  //boil::oout<<"after  "<<totalvol(clr)<<"\n";
  //std::cout.unsetf(std::ios_base::floatfield);
  //std::cout<< setprecision(6);
#ifdef DEBUG
  boil::oout<<"cipcsl2_redist:end\n";
#endif
  return;

}
