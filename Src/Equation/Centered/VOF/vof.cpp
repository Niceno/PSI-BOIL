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
  clr( *PHI.domain() ),
  nx( *PHI.domain() ),
  ny( *PHI.domain() ),
  nz( *PHI.domain() ),
  mx( *PHI.domain() ),
  my( *PHI.domain() ),
  mz( *PHI.domain() ),
  nalpha( *PHI.domain() ),
  utx( *PHI.domain() ),
  uty( *PHI.domain() ),
  utz( *PHI.domain() ),
  unliq( *PHI.domain() ),
  stmp( *PHI.domain() ),
  stmp2( *PHI.domain() ),
  stmp3( *PHI.domain() ),
  fs( *U.domain() ),
  iflag(*PHI.domain() ),
  iflagx(*PHI.domain() ),
  adens(*PHI.domain() ),
  adensgeom(*PHI.domain() ),
  heavi(&phi, NULL, &adens),
  topo(&mx,&my,&mz,&adens,&fs)

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  //kappa     = phi.shape();
  clr       = phi.shape();
  stmp      = phi.shape();
  stmp2     = phi.shape();
  stmp3     = phi.shape();
  iflag     = phi.shape();
  iflagx    = phi.shape();

  nx        = phi.shape();
  ny        = phi.shape();
  nz        = phi.shape();
  mx        = phi.shape();
  my        = phi.shape();
  mz        = phi.shape();
  nalpha    = phi.shape();
  adens     = phi.shape();
  adensgeom = phi.shape();

  utx = phi.shape();
  uty = phi.shape();
  utz = phi.shape();
  unliq = phi.shape();

  for( int b=0; b<phi.bc().count(); b++ ) {
    if(    phi.bc().type(b) == BndType::dirichlet()
        || phi.bc().type(b) == BndType::inlet()
        || phi.bc().type(b) == BndType::insert()
        || phi.bc().type(b) == BndType::convective()
       ) {
       nx.bc().type(b) = BndType::neumann();
       ny.bc().type(b) = BndType::neumann();
       nx.bc().type(b) = BndType::neumann();
       nalpha.bc().type(b) = BndType::neumann();
       adens.bc().type(b) = BndType::neumann();
       adensgeom.bc().type(b) = BndType::neumann();
       mx.bc().type(b) = BndType::neumann();
       my.bc().type(b) = BndType::neumann();
       mx.bc().type(b) = BndType::neumann();

       utx.bc().type(b) = BndType::neumann();
       uty.bc().type(b) = BndType::neumann();
       utz.bc().type(b) = BndType::neumann();
       unliq.bc().type(b) = BndType::neumann();

       boil::oout << "Adjusting b.c.s for geometrical properties at " << b
                  << boil::endl;
    }
  }

 
  mixture = flu;
  if(mixture) {
    rhol = mixt()->rho(1);
    rhov = mixt()->rho(0);
  } else {
    rhol = 1.;
    rhov = 1.;
  }

  for_m(m)
    fs(m) = (*u)(m).shape();
  bndclr = BNDCLR;

  assert(PHI.domain() == F.domain());

/* set parameters */
  //dxmin=std::min(phi.dxc(3),std::min(phi.dyc(3),phi.dzc(3)));
  dxmin = dom->dxyz_min();
  boil::cart.min_real(&dxmin);
  ww=1.0*dxmin;

  epsnorm=1.0e-12;
  phisurf=0.5;
  nlayer=6;
  tol_wall = 0.01; /* tolerance 0.99 \approx 1.0 near walls */
  curv_method = 0;
  cangle=90.0/180.0*boil::pi;
  limit_color=false;

  discretize();

  /* apply boundary condition */
  phi.bnd_update();
  phi.exchange_all();

  /* check boundary condition */
  iminp = imaxp = jminp = jmaxp = kminp = kmaxp = false; // true for periodic
  iminc = imaxc = jminc = jmaxc = kminc = kmaxc = true;  // true for cut-stencil
  iminw = imaxw = jminw = jmaxw = kminw = kmaxw = false; // true for wall
  ifull = jfull = kfull = true; // true for not a dummy direction
  // imin
  Dir d = Dir::imin();
  if (phi.bc().type_decomp(d)) {
    iminp=true;
    iminc=false;
  } else {
    if (phi.bc().type(d,BndType::periodic())) {
      iminp=true;
      iminc=false;
    } else if (phi.bc().type(d,BndType::wall())) {
      iminw=true;
    } else if (phi.bc().type(d,BndType::pseudo())) {
      iminp=true;
      iminc=false;
      ifull = false;
    }

    if (dom->bnd_symmetry(d)) iminc=false;
  }
  // imax
  d = Dir::imax();
  if (phi.bc().type_decomp(d)) {
    imaxp=true;
    imaxc=false;
  } else {
    if (phi.bc().type(d,BndType::periodic())) {
      imaxp=true;
      imaxc=false;
    } else if (phi.bc().type(d,BndType::wall())) {
      imaxw=true;
    } else if (phi.bc().type(d,BndType::pseudo())) {
      imaxp=true;
      imaxc=false;
      ifull = false;
    }
    if (dom->bnd_symmetry(d)) imaxc=false;
  }
  // jmin
  d = Dir::jmin();
  if (phi.bc().type_decomp(d)) {
    jminp=true;
    jminc=false;
  } else {
    if (phi.bc().type(d,BndType::periodic())) {
      jminp=true;
      jminc=false;
    } else if (phi.bc().type(d,BndType::wall())) {
      jminw=true;
    } else if (phi.bc().type(d,BndType::pseudo())) {
      jminp=true;
      jminc=false;
      jfull = false;
    }
    if (dom->bnd_symmetry(d)) jminc=false;
  }
  // jmax
  d = Dir::jmax();
  if (phi.bc().type_decomp(d)) {
    jmaxp=true;
    jmaxc=false;
  } else {
    if (phi.bc().type(d,BndType::periodic())) {
      jmaxp=true;
      jmaxc=false;
    } else if (phi.bc().type(d,BndType::wall())) {
      jmaxw=true;
    } else if (phi.bc().type(d,BndType::pseudo())) {
      jmaxp=true;
      jmaxc=false;
      jfull = false;
    }
    if (dom->bnd_symmetry(d)) jmaxc=false;
  }
  // kmin
  d = Dir::kmin();
  if (phi.bc().type_decomp(d)) {
    kminp=true;
    kminc=false;
  } else {
    if (phi.bc().type(d,BndType::periodic())) {
      kminp=true;
      kminc=false;
    } else if (phi.bc().type(d,BndType::wall())) {
      kminw=true;
    } else if (phi.bc().type(d,BndType::pseudo())) {
      kminp=true;
      kminc=false;
      kfull = false;
    }
    if (dom->bnd_symmetry(d)) kminc=false;
  }
  // kmax
  d = Dir::kmax();
  if (phi.bc().type_decomp(d)) {
    kmaxp=true;
    kmaxc=false;
  } else {
    if (phi.bc().type(d,BndType::periodic())) {
      kmaxp=true;
      kmaxc=false;
    } else if (phi.bc().type(d,BndType::wall())) {
      kmaxw=true;
    } else if (phi.bc().type(d,BndType::pseudo())) {
      kmaxp=true;
      kmaxc=false;
      kfull = false;
    }
    if (dom->bnd_symmetry(d)) kmaxc=false;
  }

  boil::aout<<"curv_HF::periodic= "<<boil::cart.iam()<<" "
            <<iminp<<" "<<imaxp<<" "
            <<jminp<<" "<<jmaxp<<" "
            <<kminp<<" "<<kmaxp<<"\n";

  boil::aout<<"curv_HF::wall= "<<boil::cart.iam()<<" "
            <<iminw<<" "<<imaxw<<" "
            <<jminw<<" "<<jmaxw<<" "
            <<kminw<<" "<<kmaxw<<"\n";

  boil::aout<<"curv_HF::cut-stencil= "<<boil::cart.iam()<<" "
            <<iminc<<" "<<imaxc<<" "
            <<jminc<<" "<<jmaxc<<" "
            <<kminc<<" "<<kmaxc<<"\n";

  boil::oout<<"VOF-full: "<<ifull<<" "<<jfull<<" "<<kfull<<boil::endl;

}	

/******************************************************************************/
VOF::~VOF() {
}	

