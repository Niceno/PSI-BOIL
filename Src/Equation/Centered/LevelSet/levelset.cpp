#include "levelset.h"
#include <cmath>

/******************************************************************************/
LevelSet::LevelSet(const Scalar & PHI, 
                           const Scalar & F,
                           const Vector & U, 
                           Times & T,
                           Krylov * S) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  jelly( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &jelly, NULL, S ),
  nx( *PHI.domain() ),
  ny( *PHI.domain() ),
  nz( *PHI.domain() ),
  nmag( *PHI.domain() ),
  kappa( *PHI.domain() ),
  stmp( *PHI.domain() ),
  dflag( *PHI.domain() )

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  nx      = phi.shape();
  ny      = phi.shape();
  nz      = phi.shape();
  nmag    = phi.shape();
  kappa   = phi.shape();
  stmp    = phi.shape();
  dflag   = phi.shape();

  assert(PHI.domain() == F.domain());

  /* set constants */
  pi = acos(-1.0);
  phimin =-1.0;
  phimax = 1.0;
  phisurf = 0.5*(phimin+phimax);
  dxmin = dom->dxyz_min();
  epsnorm=1.0e-12;

  /* set initial value */
  cangle=90.0;
  nredist=-1;

  insert_bc_dist(phi);
  phi.exchange_all();
}	

/******************************************************************************/
LevelSet::~LevelSet() {
}
