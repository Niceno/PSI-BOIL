#include "custom.h"
  
namespace boil {
  /******************************************************************************/
  real convective_boundary_layer_thickness(const Matter & mat,
                                           const real deltat) {
  /***************************************************************************//**
   \brief Calculate the natural convection boundary layer thickness according to
          the correlation of Kays and Crawford (Convective Heat and Mass Transfer)
      output: BL thickness
  *******************************************************************************/
    real mu      = mat.mu()->value();
    real lambda  = mat.lambda()->value();
    real cp      = mat.cp()->value(); /* J/m3 */
    real rho     = mat.rho()->value();
    real beta    = mat.beta()->value(); /* 1/T */
    real gravity = boil::g;

    return 7.14*std::pow(mu*lambda/cp/rho/gravity/deltat/beta,1./3.);
  }
}
