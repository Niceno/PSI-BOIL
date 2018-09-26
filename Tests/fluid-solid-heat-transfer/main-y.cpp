#include "Include/psi-boil.h"

/* boundary conditions */
const real LX  = 1.0;
const real LY  = 1.0; 
const real mu  = 0.01;
const real lam = 0.01;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(0.001,LX+0.001), 
             20, Periodic::no());

  Grid1D gz( Range<real>(0.0,LY),  
             20, Periodic::yes());

  Grid1D gy( Range<real>(0.0,LX), 
             20, Periodic::no());

  /*----------------+
  |  immersed body  |
  +----------------*/
  Body wall("wall.stl");

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz, & wall);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d), solid(d);
  fluid.mu    (0.001);
  fluid.lambda(0.001);
  solid.lambda(1.00);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Times time(100, 0.075); /* ndt, dt */
	
  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar t(d),   g(d);   // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  t.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(),  1.0 ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), -1.0 ) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Vector uvw(d);                          // velocity
  Enthalpy enth(t, g, uvw, time, solver, &fluid, &solid);

//..  AC multigrid( &enth );

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    OMS(A);
    enth.new_time_step();
    OMS(AA);
    enth.solve(ResRat(0.0001), "temperature");
    OMS(AAA);
//..    multigrid.vcycle(ResRat(1e-4));
  }
  boil::plot->plot(t,   "t",  time.current_step());

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(t, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-y.cpp,v 1.3 2011/05/25 11:13:59 niceno Exp $'/
+-----------------------------------------------------------------------------*/
