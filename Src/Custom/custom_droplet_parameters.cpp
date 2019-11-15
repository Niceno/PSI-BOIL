#include "custom.h"

namespace boil {
  void droplet_parameters_2D(const real cang, const real area,
                             real & radius, real & zcent, real & chord) {
    real cangle = cang * boil::pi/180.;

    if(cangle > 0.5*boil::pi) {
      real alpha = 2.*(boil::pi-cangle);
      radius = sqrt(area/(boil::pi-0.5*(alpha-sin(alpha))));
      zcent  = radius*cos(0.5*alpha);
      chord  = 2.*radius*sin(0.5*alpha);
    } else {
      real alpha = 2.*cangle;
      radius = sqrt(area/(0.5*(alpha-sin(alpha))));
      zcent  = -radius*cos(0.5*alpha);
      chord  = 2.*radius*sin(0.5*alpha);
    }

    return;
  }

  void droplet_parameters_3D(const real cang, const real volume,
                             real & radius, real & zcent, real & chord) {
    real cangle = cang * boil::pi/180.;

    if(cangle > 0.5*boil::pi) {
      real alpha = 2.*(boil::pi-cangle);
      real nu = 1. - cos(alpha/2.);
      radius = pow(volume/(boil::pi/3.) / (4. + nu*nu*nu - 3. * nu*nu), 1./3.);
      zcent  = radius*cos(0.5*alpha);
      chord  = 2.*radius*sin(0.5*alpha);
    } else {
      real alpha = 2.*cangle;
      real nu = 1. - cos(alpha/2.);
      radius = pow(volume/(boil::pi/3.) / (3.* nu*nu - nu*nu*nu), 1./3.);
      zcent  = -radius*cos(0.5*alpha);
      chord  = 2.*radius*sin(0.5*alpha);
    }

    return;
  }
}
