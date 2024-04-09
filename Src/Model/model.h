#ifndef MODEL_H
#define MODEL_H

#include <cassert>

#include "../Equation/Staggered/Momentum/momentum.h"

/////////////
//         //
//  Model  //
//         //
/////////////
class Model {

  public:
    void smagorinsky( const Momentum * m, Scalar * mu_t, real c_s,
                      const real * dist_max=NULL,
                      const Scalar * dist = NULL, Scalar * yplus = NULL) const;
    void wale       ( const Momentum * m, Scalar * mu_t, real c_w ) const;
    void smagorinsky_twophase( const Momentum * m, Scalar * mu_t, real c_s,
                      const Scalar * c,
                      const real * dist_max=NULL,
                      const Scalar * dist = NULL, Scalar * yplus = NULL) const;
    void tau_wall   ( Momentum * m, const Scalar & dist,
                      Vector * fbnd, Scalar * y_plus = NULL) const;

  protected:
    real tangential_velocity( const real u, const real v, const real w,
                              const real nx, const real ny, const real nz,
                              real * tx, real * ty, real * tz) const;
    void wall_function( const Vector & uvw, const Matter & fluid,
                        const Scalar & dist, 
                        real * tau_w, real * y_pl,
                        const int i, const int j, const int k,
                        real * tx, real * ty, real * tz ) const;
};

#endif
