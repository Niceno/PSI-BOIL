#include "vof.h"
#include <cmath>

/******************************************************************************/
VOF::VOF(const Scalar & PHI, 
         const Scalar & F,
         const Scalar & K,
         const Vector & U, 
         Times & T,
         Krylov * S,
         Vector * BNDCLR,
         Matter * flu) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  jelly( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &jelly, NULL, S ),
  kappa( &K ),
  nx( *PHI.domain() ),
  ny( *PHI.domain() ),
  nz( *PHI.domain() ),
  mx( *PHI.domain() ),
  my( *PHI.domain() ),
  mz( *PHI.domain() ),
  unliq( *PHI.domain() ),
  utliq( *PHI.domain() ),
  utx( *PHI.domain() ),
  uty( *PHI.domain() ),
  utz( *PHI.domain() ),
  uliq ( *U   .domain() ),
  nalpha( *PHI.domain() ),
  nmag( *PHI.domain() ),
  stmp( *PHI.domain() ),
  stmp2( *PHI.domain() ),
  stmp3( *PHI.domain() ),
  fs( *U.domain() ),
  fsx( *PHI.domain() ),
  fsy( *PHI.domain() ),
  fsz( *PHI.domain() ),
  iflag(*PHI.domain() ),
  iflagx(*PHI.domain() ),
  iflagy(*PHI.domain() ),
  iflagz(*PHI.domain() ),
  adens(*PHI.domain() ),
  adensgeom(*PHI.domain() )

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  nx        = phi.shape();
  ny        = phi.shape();
  nz        = phi.shape();
  mx        = phi.shape();
  my        = phi.shape();
  mz        = phi.shape();
  unliq     = phi.shape();
  utliq     = phi.shape();
  for_m(m) uliq(m) = U(m).shape();
  utx     = phi.shape();
  uty     = phi.shape();
  utz     = phi.shape();
  nmag      = phi.shape();
  nalpha    = phi.shape();
  stmp      = phi.shape();
  stmp2     = phi.shape();
  stmp3     = phi.shape();
  fsx       = phi.shape();
  fsy       = phi.shape();
  fsz       = phi.shape();
  iflag     = phi.shape();
  iflagx    = phi.shape();
  iflagy    = phi.shape();
  iflagz    = phi.shape();
 
  mixture = flu;
  if(mixture) {
    rhol = mixt()->rho(1);
    rhov = mixt()->rho(0);
  } else {
    rhol = 1.;
    rhov = 1.;
  }

  adens = phi.shape();
  adensgeom = phi.shape();
  for_m(m) {
    fs(m) = (*u)(m).shape();
  }

  bndclr = BNDCLR;

  assert(PHI.domain() == F.domain());

  pi = acos(-1.0);
  //dxmin=std::min(phi.dxc(3),std::min(phi.dyc(3),phi.dzc(3)));
  dxmin = dom->dxyz_min();
  boil::cart.min_real(&dxmin);
  ww=1.0*dxmin;

  epsnorm=1.0e-12;
  phisurf=0.5;
  tol_wall = 0.01; /* tolerance 0.99 \approx 1.0 near walls */
  tol_flux = 3e-3; /* during inner iterations */
  flux_cfl = 0.2;  /* used in case 3 flux calculations */
  maxiter = 10;    /* maximal number of iterations */

  /* set initial value */
  nlayer=4;
  n_ext_fs=5;

  //alloc3d(& iflag, phi.ni(), phi.nj(), phi.nk());

#if 0
  f_w=0.0;
  f_e=0.0;
  f_t=0.0;
  f_b=0.0;
  f_n=0.0;
  f_s=0.0;
#endif

  discretize();

  /* apply boundary condition */
  phi.bnd_update();
  phi.exchange_all();

}	

/******************************************************************************/
VOF::~VOF() {
}	

