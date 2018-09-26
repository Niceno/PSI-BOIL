#include "model.h"

/******************************************************************************/
real Model::tangential_velocity( const real ux, const real uy, const real uz,
                                 const real nx, const real ny, const real nz, 
                                 real     * tx, real     * ty, real     * tz) 
                                 const {
/*-----------------------------------------------------+ 
|  tangential vector: t = n x (u x n) / |n x (u x n)|  |
+-----------------------------------------------------*/

  /* q = u x n */
  const real qx =    uy*nz - uz*ny;
  const real qy = - (ux*nz - uz*nx);
  const real qz =    ux*ny - uy*nx;

  /* t = n x q = n x (u x n) */
  * tx =    ny*qz - nz*qy;
  * ty = - (nx*qz - nz*qx);
  * tz =    nx*qy - ny*qx;

  /* tm = | u x (u x n) | */
  real tm = sqrt( (*tx)*(*tx) + (*ty)*(*ty) + (*tz)*(*tz) );
  tm = std::max(tm,boil::atto);  // avoid divided-by-zero

  * tx /= tm;
  * ty /= tm;
  * tz /= tm;

  return *tx*ux + *ty*uy + *tz*uz;

  /*-----------------------------+
  |  old: simple and stupid      |
  real un = u*nx + v*ny + w*nz;
  real ua2 = u*u + v*v + w*w;

  return sqrt(ua2 - un);
  +-----------------------------*/
}

/*-----------------------------------------------------------------------------+
 '$Id: model_tangential_velocity.cpp,v 1.3 2014/10/28 12:51:52 sato Exp $'/
+-----------------------------------------------------------------------------*/
