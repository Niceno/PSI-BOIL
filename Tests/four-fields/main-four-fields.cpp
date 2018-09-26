/*----------------------+
|                       |
|  set lcnd to 0.5 ...  |
|                       |
+----------------------*/
#include "Include/psi-boil.h"

/* boundary conditions */
const int NX = 96;
const int NZ=  NX * 3/2;

const real D0 = 0.0005;
const real D1 = 0.00075;
const real D2 = 0.001;
const real D3 = 0.00125;

const real LX = D3*8.0;           
const real LZ = LX*NZ/NX; 

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC(AsNodes::yes());

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), NX, Periodic::no());
  Grid1D gz( Range<real>(-LZ/2.0, LZ/2.0), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  Times time(2560, 0.00004); /* ndt, dt */
	
  const int         nsint = 100;
        int         csint =   0;
  std::vector<real> save_instants;
  save_instants.resize(nsint+1); /* to store last too */
  for(int i=0; i<=nsint; i++) {
    save_instants[i] = (real)i * time.total_time() / (real)nsint;
    OPR(save_instants[i]);
  }

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::ic3());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p    (d, "press"),     f  (d); // p.
  Scalar c_sep(d, "separated"), g  (d); // separated phase
  Scalar c_bub(d, "bubbles");           // dispersed phase
  Scalar c_dro(d, "droplets");          // dispersed phase
  Scalar press(d, "pressure");

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  press = p.shape();
  c_sep = p.shape();
  c_bub = p.shape();
  c_dro = p.shape();
  
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (   0.00001);
  air  .rho   (   1.2    );
  water.mu    (   0.001  );
  water.rho   ( 999.0    );

  Matter mixed(water, air, c_sep);
  mixed.sigma (   0.07   );

  ColorFunction conc(c_sep,   g, 1.0, 1.0, uvw, time, solver); 
  conc.convection_set(ConvScheme::superbee());
  conc.sharpen();

  Dispersed bubbles (c_bub, & c_sep, 0, uvw, time, &mixed); 
  Dispersed droplets(c_dro, & c_sep, 1, uvw, time, &mixed); 

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());

  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);

  /*====================+ 
  |                     |
  |  initial condition  |
  |                     |
  +====================*/

  /*-------------+ 
  |  all liquid  |
  +-------------*/
  c_sep = 1.0;

  /*---------+
  |  circle  |
  +---------*/
//for_vijk(c_sep,i,j,k) {
//  real dist = sqrt( c_sep.xc(i)*c.xc(i) + c.yc(j)*c.yc(j) + c.zc(k)*c.zc(k) ); 
//  if(dist < RB)
//    c_sep[i][j][k] = 0.0; /* gas   */
//  else
//    c_sep[i][j][k] = 1.0; /* liquid*/
//}

  /*------------+
  |  interface  |
  +------------*/
  for_vijk(c_sep,i,j,k) {
    if(c_sep.zc(k) > LZ/8)
      c_sep[i][j][k] = 0.0;   /* gas   */
  }

  c_sep.exchange();

  /*-------------------+
  |  dispsersed phase  |
  +-------------------*/
  const real offst = 2.0 * D3;
  bubbles.add( Particle( Position( -offst*.707, -offst*.707, -LZ/2.0 + 3.5*D3), Diameter(D0) ) );
  bubbles.add( Particle( Position(  offst*.707, -offst*.707, -LZ/2.0 + 3.0*D3), Diameter(D1) ) );
  bubbles.add( Particle( Position(  offst*.707,  offst*.707, -LZ/2.0 + 2.5*D3), Diameter(D2) ) );
  bubbles.add( Particle( Position( -offst*.707,  offst*.707, -LZ/2.0 + 2.0*D3), Diameter(D3) ) );
  c_bub.exchange();
  droplets.add( Particle( Position(  offst,    0.0,  LZ/2.0 - 3.5*D3), Diameter(D0) ) );
  droplets.add( Particle( Position(    0.0,  offst,  LZ/2.0 - 3.0*D3), Diameter(D1) ) );
  droplets.add( Particle( Position( -offst,    0.0,  LZ/2.0 - 2.5*D3), Diameter(D2) ) );
  droplets.add( Particle( Position(    0.0, -offst,  LZ/2.0 - 2.0*D3), Diameter(D3) ) );
  c_dro.exchange();

  boil::plot->plot(uvw, c_sep, c_bub, c_dro, "uvw-phases", 0);
//boil::plot->plot(*(mixed.rho()), "density",  time.current_step());

  /*============+ 
  |             |
  |  time loop  |
  |             |
  +============*/
  for(time.start(); time.end(); time.increase()) {

    /* advance */
    conc.new_time_step();
    conc.convection();
    conc.advance();
    conc.sharpen();
    conc.tension(&xyz, mixed);
//  boil::plot->plot(xyz,  c_sep, "xyz,c",  time.current_step());
   
    const Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] -= 9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    bubbles .advance(& xyz);
    droplets.advance(& xyz);

    ns.cfl_max();

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();

    ns.grad(press);
    ns.convection();
    ns.solve(ResRat(0.01));

    p = 0.0;
    p.exchange();
    multigrid.vcycle(ResRat(1e-2));
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();

    /* time step control */
    time.control_dt(ns.cfl_max(), 0.33, 1.0);

    /*-------+
    |  plot  |
    +-------*/
    if(csint <= nsint) {
      if(time.current_time() > save_instants[csint]) {
        boil::plot->plot(uvw, c_sep, c_bub, c_dro, "uvw-phases",  csint);
        csint++;
      }
    }

    /*---------+
    |  exit ?  |
    +---------*/
    std::ifstream infile;
    infile.open ("stop.now", std::ifstream::in);
    if(infile.good()) exit(0);
    infile.close();
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
//boil::plot->plot(uvw,c_sep, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-four-fields.cpp,v 1.2 2013/08/19 07:20:51 niceno Exp $'/
+-----------------------------------------------------------------------------*/
