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
  stmp4( *PHI.domain() ),
  stmp5( *PHI.domain() ),
  stmp6( *PHI.domain() ),
  fs( *U.domain() ),
  fluxmax( *U.domain() ),
  sosflux( *U.domain() ),
  iflag(*PHI.domain() ),
  iflagx(*PHI.domain() ),
  iflagy(*PHI.domain() ),
  iflagz(*PHI.domain() ),
  adens(*PHI.domain() ),
  //adensgeom(*PHI.domain() ),
  heavi(&phi, NULL, &adens),
  topo(&mx,&my,&mz,&adens,&fs)

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
  stmp4     = phi.shape();
  stmp5     = phi.shape();
  stmp6     = phi.shape();
  iflag     = phi.shape();
  iflagx    = phi.shape();
  iflagy    = phi.shape();
  iflagz    = phi.shape();
  adens     = phi.shape();
 
  mixture = flu;
  if(mixture) {
    rhol = mixt()->rho(1);
    rhov = mixt()->rho(0);
  } else {
    rhol = 1.;
    rhov = 1.;
  }

  //adensgeom = phi.shape();
  for_m(m) {
    fs(m) = (*u)(m).shape();
    sosflux(m) = (*u)(m).shape();
    fluxmax(m) = (*u)(m).shape();
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
  nlayer=6;
  n_ext_fs=5;
  tol_wall = 0.01; /* tolerance 0.99 \approx 1.0 near walls */
  tol_flux = 3e-3; /* during inner iterations */
  tol_ext = 1e-7; /* extrapolation tolerance */
  flux_cfl = 0.2;  /* used in case 3 flux calculations */
  maxiter = 10;    /* maximal number of iterations */
  curv_method = 0;

  discretize();

  /* apply boundary condition */
  phi.bnd_update();
  phi.exchange_all();

  /* check boundary condition */
  iminp = imaxp = jminp = jmaxp = kminp = kmaxp = false;
  // imin
  Dir d = Dir::imin();
  if (phi.bc().type_decomp(d)) iminp=true;
  if (phi.bc().type(d,BndType::periodic())) iminp=true;
  // imax
  d = Dir::imax();
  if (phi.bc().type_decomp(d)) imaxp=true;
  if (phi.bc().type(d,BndType::periodic())) imaxp=true;
  // jmin
  d = Dir::jmin();
  if (phi.bc().type_decomp(d)) jminp=true;
  if (phi.bc().type(d,BndType::periodic())) jminp=true;
  // jmax
  d = Dir::jmax();
  if (phi.bc().type_decomp(d)) jmaxp=true;
  if (phi.bc().type(d,BndType::periodic())) jmaxp=true;
  // kmin
  d = Dir::kmax();
  if (phi.bc().type_decomp(d)) kminp=true;
  if (phi.bc().type(d,BndType::periodic())) kminp=true;
  // kmax
  d = Dir::kmax();
  if (phi.bc().type_decomp(d)) kmaxp=true;
  if (phi.bc().type(d,BndType::periodic())) kmaxp=true;

  boil::aout<<"curv_HF::Bnd= "<<iminp<<" "<<imaxp<<" "<<jminp<<" "<<jmaxp<<" "
            <<kminp<<" "<<kmaxp<<"\n";

}	

/******************************************************************************/
VOF::~VOF() {
}	

