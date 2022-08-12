#include "cipcsl2.h"
#include <cmath>
//#define DEBUG

/******************************************************************************/
CIPCSL2::CIPCSL2(const Scalar & PHI, 
                 const Scalar & F,
                 const Scalar & K,
                 const Vector & U, 
                 Times & T,
                 Linear * S) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  jelly( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &jelly, NULL, S ),
  clr   ( *PHI.domain() ),
  sclr  ( *PHI.domain() ),
  nx    ( *PHI.domain() ),
  ny    ( *PHI.domain() ),
  nz    ( *PHI.domain() ),
  atmp  ( *PHI.domain() ),
  stmp  ( *PHI.domain() ),
  iflag ( *PHI.domain() ),
  wflag ( *PHI.domain() ),
  intflag( *PHI.domain() ),
  scheme( PHI ),
  fn    ( *PHI.domain() ),
  dist  ( *PHI.domain() ),
  alp   ( *PHI.domain() ),
  sxyz  ( *U.domain() ),
  fs    ( *U.domain() ),
  adens ( *PHI.domain() ),
  kappa ( &K )
  //heavi(&phi, NULL, &adens),
  //topo(&phi,&phi,&nx,&ny,&nz,&adens,&fs,&intflag)

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  clr     = phi.shape();
  sclr    = phi.shape();
  nx      = phi.shape();
  ny      = phi.shape();
  nz      = phi.shape();
  atmp    = phi.shape();
  stmp    = phi.shape();
  iflag   = phi.shape();
  wflag   = phi.shape();
  intflag = phi.shape();
  fn      = phi.shape();
  dist    = phi.shape();
  alp     = phi.shape();
  adens   = phi.shape();
  kappa   = phi.shape();

  for_m(m){
    sxyz(m) = (*u)(m).shape();
    fs(m)   = (*u)(m).shape();
  }
  assert(PHI.domain() == F.domain());

  /* set constants */
  phimin = 0.0;
  phimax = 1.0;
  phisurf = 0.5*(phimin+phimax);
  dxmin = dom->dxyz_min();
  ww=1.0*dxmin; // default value for ww
  tol_wall = 0.01;

  /* runtime polymorphism */
  heavi = new MarchingCubes(&phi,NULL,&adens);
  topo = new Topology(&phi,&phi,&nx,&ny,&nz,&adens,&fs,&intflag,phisurf);

  epsnorm=1.0e-12;
  eps_clr = 1.0e-4; // epsilon for color function

  /* set initial value */
  nredist=1;
  itsharpen=4;
  ialpcal=0;
  cangle=90.0/180.0*boil::pi;
  itsmear=10;
  eps_st=1.5;
  sum_outlet=0.0;
  sum_outletm=0.0;
  nlayer=16;
  localSharpen=true;
  use_dist_for_kappa=true;
  extrapolate_ib=true;

  /* allocate array */
  alloc3d(& vel,    phi.ni()+1, phi.nj()+1, phi.nk()+1);
  alloc3d(& delrho, phi.ni()+1, phi.nj()+1, phi.nk()+1);

#ifdef DEBUG
  std::cout<<"cipcsl2:alloc \n";
#endif

  //discretize();
#ifdef DEBUG
  std::cout<<"cipcsl2:discretize \n";
#endif

  /* set in init() */
  is_initialized = false;

  init();
#ifdef DEBUG
  std::cout<<"cipcsl2:init \n";
#endif
}	

/******************************************************************************/
CIPCSL2::~CIPCSL2() {
}
