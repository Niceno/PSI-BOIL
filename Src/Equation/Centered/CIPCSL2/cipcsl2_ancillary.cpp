#include "cipcsl2.h"
//#define DEBUG

/******************************************************************************/
void CIPCSL2::ancillary() {
/***************************************************************************//**
*  \brief Calculate ancillary vof parameters.
*******************************************************************************/

  boil::timer.start("cipcsl2 ancillary");

  /* wall values extrapolation */
  update_at_walls(phi);
#ifdef DEBUG
  std::cout<<"ancillary:update_at_walls\n";
#endif

  /* calculate free surface position,
     no subgrid interfaces near walls! */
  heavi->cal_fs_interp(phi,fs,0.0,false);
#ifdef DEBUG
  std::cout<<"ancillary:heavi\n";
#endif

  /* normal vector */
  distfunc(phi,24); 
#ifdef DEBUG
  std::cout<<"ancillary:distfunc\n";
#endif
  gradphic(dist);
#ifdef DEBUG
  std::cout<<"ancillary:gradphic\n";
#endif

  /* calculate area density */
  heavi->calculate_adens();

  /* flag the interface */
  interfacial_flagging(phi);

  boil::timer.stop("cipcsl2 ancillary");

  return;
}

