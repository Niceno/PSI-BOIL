#include "solidparticle.h"

/***************************************************************************//**
*  Constructor for solidparticle
*******************************************************************************/
SolidParticle::SolidParticle ( const Vector & v, const Times * t,
                               Matter * mat,
                               const real gx, const real gy, const real gz,
                               const Scalar * sca1,
                               const Scalar * sca2,
                               const Scalar * sca3) :
  Pathline(v, t, sca1, sca2, sca3),
  flu(mat) 
{
  gravity_x = gx;
  gravity_y = gy;
  gravity_z = gz;
  b_dia_den = true;
}

/******************************************************************************/
//SolidParticle::~SolidParticle() {
//}
