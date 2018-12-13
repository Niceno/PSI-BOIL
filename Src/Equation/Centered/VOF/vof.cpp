#include "vof.h"
#include <cmath>

/******************************************************************************/
VOF::VOF(const Scalar & PHI, 
         const Scalar & F,
         const Scalar & K,
         const Vector & U, 
         Times & T,
         Krylov * S,
         Vector * BNDCLR) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  jelly( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &jelly, NULL, S ),
  kappa( &K ),
  clr( *PHI.domain() ),
  nx( *PHI.domain() ),
  ny( *PHI.domain() ),
  nz( *PHI.domain() ),
  nalpha( *PHI.domain() ),
  nmag( *PHI.domain() ),
  clrn( *PHI.domain() ),
  stmp( *PHI.domain() ),
  fs( *U.domain() ),
  fsx( *PHI.domain() ),
  fsy( *PHI.domain() ),
  fsz( *PHI.domain() ),
  iflag(*PHI.domain() ),
  iflagx(*PHI.domain() ),
  iflagy(*PHI.domain() ),
  iflagz(*PHI.domain() ),
  adens(*PHI.domain() )

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  clr       = phi.shape();
  nx        = phi.shape();
  ny        = phi.shape();
  nz        = phi.shape();
  nmag      = phi.shape();
  nalpha    = phi.shape();
  stmp      = phi.shape();
  fsx       = phi.shape();
  fsy       = phi.shape();
  fsz       = phi.shape();
  iflag     = phi.shape();
  iflagx    = phi.shape();
  iflagy    = phi.shape();
  iflagz    = phi.shape();
 
  adens = phi.shape();
  for_m(m) {
    fs(m) = (*u)(m).shape();
  }

  bndclr = BNDCLR;

  assert(PHI.domain() == F.domain());

  pi = acos(-1.0);
  dxmin=std::min(phi.dxc(1),std::min(phi.dyc(1),phi.dzc(1)));
  boil::cart.min_real(&dxmin);
  ww=1.0*dxmin;

  epsnorm=1.0e-12;
  phisurf=0.5;

/* set initial value */
  nlayer=4;
  n_ext_fs=5;

  //alloc3d(& iflag, phi.ni(), phi.nj(), phi.nk());

  discretize();

  /* apply boundary condition */
  phi.bnd_update();
  phi.exchange_all();

}	

/******************************************************************************/
VOF::~VOF() {
}	

