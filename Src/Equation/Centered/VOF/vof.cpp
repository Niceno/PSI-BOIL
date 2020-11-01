#include "vof.h"
#include <cmath>

/******************************************************************************/
VOF::VOF(const Scalar & PHI, 
         const Scalar & F,
         const Scalar & K,
         const Vector & U, 
         Times & T,
         Linear * S,
         Vector * BNDCLR) :
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
  nalpha( *PHI.domain() ),
  stmp( *PHI.domain() ),
  stmp2(*PHI.domain() ),
  pold_neg(*PHI.domain() ),
  pold_pos(*PHI.domain() ),
  fs( *U.domain() ),
  vflow( *U.domain() ),
  iflag(*PHI.domain() ),
  tempflag(*PHI.domain() ),
  tempflag2(*PHI.domain() ),
  adens(*PHI.domain() ),
  //topo(&phi,&color(),&mx,&my,&mz,&adens,&fs,&iflag),
  norm_method_advance(NormMethod::Mixed()),
  norm_method_curvature(NormMethod::Young()),
  mcomp_for_elvira(Comp::undefined()), /* undefined for 3D */
  bulk_curv_method(CurvMethod::HF()),
  wall_curv_method(CurvMethod::DivNorm()),
  subgrid_method(SubgridMethod::PLIC()),
  topo_method(TopoMethod::Hybrid()),
  hf_set(),
  bflag_struct(PHI)

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
#if 0 /* don't use this, it creates BndCnd pointers */
  nx     = phi.shape();
  ny     = phi.shape();
  nz     = phi.shape();
  mx     = phi.shape();
  my     = phi.shape();
  mz     = phi.shape();
  nalpha = phi.shape();
  stmp   = phi.shape();
  stmp2  = phi.shape();
  iflag  = phi.shape();
  adens  = phi.shape();
  pold_neg  = phi.shape();
  pold_pos  = phi.shape();
  tempflag  = phi.shape();
  tempflag2 = phi.shape();
#else
  nx.copy_shape(phi.shape());
  ny.copy_shape(phi.shape());
  nz.copy_shape(phi.shape());
  mx.copy_shape(phi.shape());
  my.copy_shape(phi.shape());
  mz.copy_shape(phi.shape());
  nalpha.copy_shape(phi.shape());
  stmp  .copy_shape(phi.shape());
  stmp2 .copy_shape(phi.shape());
  iflag .copy_shape(phi.shape());
  adens .copy_shape(phi.shape());
  pold_neg .copy_shape(phi.shape());
  pold_pos .copy_shape(phi.shape());
  tempflag .copy_shape(phi.shape());
  tempflag2.copy_shape(phi.shape());
#endif
  
  phisurf=0.5;
  topo = new Topology(&phi,&color(),&mx,&my,&mz,&adens,&fs,&iflag,phisurf);

  /* dangerous: needs to be overridden in derived class!!! */
  topo->clr = &color();

  for( int b=0; b<phi.bc().count(); b++ ) {
    if(    phi.bc().type(b) == BndType::dirichlet()
        || phi.bc().type(b) == BndType::inlet()
        || phi.bc().type(b) == BndType::outlet()
        || phi.bc().type(b) == BndType::insert()
        || phi.bc().type(b) == BndType::convective()
       ) {
       nx.bc().type(b) = BndType::neumann();
       ny.bc().type(b) = BndType::neumann();
       nx.bc().type(b) = BndType::neumann();
       nalpha.bc().type(b) = BndType::neumann();
       adens.bc().type(b) = BndType::neumann();
       mx.bc().type(b) = BndType::neumann();
       my.bc().type(b) = BndType::neumann();
       mx.bc().type(b) = BndType::neumann();
       iflag.bc().type(b) = BndType::neumann();
       tempflag.bc().type(b) = BndType::neumann();
       tempflag2.bc().type(b) = BndType::neumann();
       stmp.bc().type(b) = BndType::neumann();
       stmp2.bc().type(b) = BndType::neumann();
       pold_neg.bc().type(b) = BndType::neumann();
       pold_pos.bc().type(b) = BndType::neumann();

       boil::oout << "Adjusting b.c.s for geometrical properties at " << b
                  << boil::endl;
    }
  }

  /* used in ev_solve */
  pold_neg = 0.;
  pold_pos = 0.;

  for_m(m) {
    fs(m) = (*u)(m).shape();
    vflow(m) = (*u)(m).shape();
  }

  for_m(m)
    for_avmijk(fs,m,i,j,k)
      fs[m][i][j][k] = boil::unreal;

  bndclr = BNDCLR;

  assert(PHI.domain() == F.domain());

  /* runtime polymorphism */
  if(phi.domain()->is_axisymmetric()) {
    /* dangerous: needs to be overridden in the derived class!!! */
    heavi = new MSaxisym(&color(),NULL,&adens);
  } else {
    if     (phi.bc().type(Dir::imin(),BndType::pseudo()))
      heavi = new MarchingSquares(Comp::i(),&color(),NULL,&adens);
    else if(phi.bc().type(Dir::jmin(),BndType::pseudo()))
      heavi = new MarchingSquares(Comp::j(),&color(),NULL,&adens);
    else if(phi.bc().type(Dir::kmin(),BndType::pseudo()))
      heavi = new MarchingSquares(Comp::k(),&color(),NULL,&adens);
    else
      heavi = new MarchingCubes(&color(),NULL,&adens);
  }

  /* set in init() */
  is_initialized = false;

  /* set parameters */
  //dxmin=std::min(phi.dxc(3),std::min(phi.dyc(3),phi.dzc(3)));
  dxmin = dom->dxyz_min();
  boil::cart.min_real(&dxmin);
  ww=1.0*dxmin;

  epsnorm=1.0e-12;
  tol_wall = 0.01; /* tolerance 0.99 \approx 1.0 near walls */
  cangle=90.0/180.0*boil::pi;
  mult_wall = 1.;
  Nfilm_crit = boil::unint;
  limit_color=false;
  use_interp=false;
  store_pressure_extrap=false;
  niter_pressure_extrap=1000;

  discretize();

  /* apply boundary condition */
  phi.bnd_update();
  phi.exchange_all();

}	

/******************************************************************************/
VOF::~VOF() {
  delete heavi;
}	

