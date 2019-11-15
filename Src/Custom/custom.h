#ifndef CUSTOM_H
#define CUSTOM_H

#include "../Field/Scalar/scalar.h"
#include <vector>

////////////////////////
//                    //
//  Custom functions  //
//                    //
////////////////////////
namespace boil {
  /* droplet parameters */
  void droplet_parameters_2D(const real cang, const real area,
                             real & radius, real & zcent, real & chord);
  void droplet_parameters_3D(const real cang, const real volume,
                             real & radius, real & zcent, real & chord);

  /* scalar error */
  real l2_scalar_error(const Scalar & sca, const Scalar & scb);
  real li_scalar_error(const Scalar & sca, const Scalar & scb);

  /* setup circle */
  void setup_circle_xz(Scalar & c, const real radius, const real xcent, const real zcent);
}

#endif
