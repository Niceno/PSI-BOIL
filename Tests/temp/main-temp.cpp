#include "Include/psi-boil.h"

using namespace std;

/* boundary conditions */
const real LX  = 1.0;
const real LY  = 0.25;
const real mu  = 0.01;
const real lam = 0.01;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(0,LX), 
             Range<real>( LX/256.0, LX/256.0 ),
             128, Periodic::no());

  Grid1D gy( Range<real>(0,LY),  
             Range<real>( LY/16.0, LY/16.0 ),
             16, Periodic::yes());

  Grid1D gz( Range<real>(-LX/2.0, LX/2.0), 
             Range<real>( LX/256.0, LX/256.0 ),
             128, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter liquid(d);
  liquid.mu    (0.01);
  liquid.lambda(0.01);
  liquid.rho   (1.00);

  liquid.variable( Set::mu    () );
  liquid.variable( Set::lambda() );
  liquid.variable( Set::rho   () );
  liquid.variable( Set::cp    () );

  Matter solid (d);
  solid.lambda(1.00);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Times time(400, 0.020); /* ndt, dt */

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p(d),   f(d);   // pressure
  Scalar t(d),   g(d);   // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  
  t.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::imin(), BndType::convective(), -1.0, 0.2) );
  t.bc().add( BndCnd( Dir::imax(), BndType::convective(),  1.0, 0.2) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns  (uvw, xyz,      time, solver, & liquid);
  Enthalpy enth(t,   g,   uvw, time, solver, & liquid);
  Pressure pr  (p,   f,   uvw, time, solver, & liquid);

  AC multigrid( &pr );

//enth.convection_set(TimeScheme::backward_euler());
//enth.diffusion_set (TimeScheme::backward_euler());

//uvw.load("uvw");

  for(time.start(); time.end(); time.increase()) {

    /* energy */
    enth.discretize();
    enth.new_time_step();
    enth.convection();
    enth.solve(ResRat(0.0001), "temperature");

    /* momentum */
    const Comp m = Comp::w();
    for(int i=1; i<t.ni()-1; i++)
      for(int j=1; j<t.nj()-1; j++)
        for(int k=2; k<t.nk()-1; k++)
          xyz[m][i][j][k] = t[i][j][k]*xyz.dV(m,i,j,k)*liquid.rho(m,i,j,k);

    ns.discretize();
    ns.cfl_max();
    ns.new_time_step();
    ns.solve(ResRat(0.0001));

    /* pressure */
    for(int i=0; i<p.ni(); i++)
      for(int j=0; j<p.nj(); j++)
        for(int k=0; k<p.nk(); k++)
          p[i][j][k] = 0.0;
    
    multigrid.vcycle(ResRat(0.01));
    p.exchange();

    /* momentum correction */
    ns.project(p);
    uvw.exchange();

    if(time.current_step() % 25 == 0) {
      boil::plot->plot(uvw,"uvw", time.current_step());
      boil::plot->plot(p,   "p",  time.current_step());
      boil::plot->plot(t,   "t",  time.current_step());
    }

  }
  boil::plot->plot(uvw,"uvw", time.current_step());
  boil::plot->plot(p,   "p",  time.current_step());
  boil::plot->plot(t,   "t",  time.current_step());

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(t, "test", 0);
}
