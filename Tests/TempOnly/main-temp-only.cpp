#include "Include/psi-boil.h"

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
             Range<real>( LX/128.0, LX/128.0 ),
             128, Periodic::no());

  Grid1D gy( Range<real>(0,LY),  
             Range<real>( LY/16.0, LY/16.0 ),
             16, Periodic::yes());

  Grid1D gz( Range<real>(0,LX), 
             Range<real>( LX/128.0, LX/128.0 ),
             128, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu    (0.01);
  fluid.lambda(0.01);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Times time(10, 0.075); /* ndt, dt */
	
  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar t(d),   g(d);   // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  t.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );

  t.bc().add( BndCnd( Range<int>(  1, 64),
                      Range<int>(t.sj(),t.ej()),
                      Dir::kmin(), BndType::dirichlet(), -1.0));
  t.bc().add( BndCnd( Range<int>( 65,128),
                      Range<int>(t.sj(),t.ej()),
                      Dir::kmin(), BndType::dirichlet(), +1.0));
  t.bc().add( BndCnd( Range<int>(  1, 64),
                      Range<int>(t.sj(),t.ej()),
                      Dir::kmax(), BndType::dirichlet(), +1.0));
  t.bc().add( BndCnd( Range<int>( 65,128),
                      Range<int>(t.sj(),t.ej()),
                      Dir::kmax(), BndType::dirichlet(), -1.0));

  t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Vector uvw(d);                          // velocity
  Enthalpy enth(t, g, uvw, time, solver, &fluid);

  AC multigrid( &enth );
//uvw.load("uvw");

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    enth.new_time_step();
    // enth.solve(ResRat(0.0001), "temperature");
    multigrid.vcycle(ResRat(1e-4));

    if(time.current_step() % 60 == 0) {
      boil::plot->plot(t,   "t",  time.current_step());
    }

  }
  boil::plot->plot(t,   "t",  time.current_step());

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(t, "test", 0);
}
