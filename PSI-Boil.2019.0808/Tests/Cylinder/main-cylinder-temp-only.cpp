#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*----------------+
  |  immersed body  |
  +----------------*/
  Body cyl("cylinder_128.stl");
  boil::plot->plot(cyl, "cylinder");

  /*-----------------------+
  |  computational domain  |
  +-----------------------*/
  Grid1D gx( Range<real>(-1.0001, 1.0001),  64, Periodic::no() ); 
  Grid1D gy( Range<real>(-1.0001, 1.0001),  64, Periodic::yes() ); 
  Grid1D gz( Range<real>(-1.0001, 1.0001),  64, Periodic::no() ); 

  Domain d(gx, gy, gz, & cyl);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d), solid(d);
  fluid.mu    (0.01);
  fluid.lambda(0.01);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Times time(20, 20.0); /* ndt, dt */
	
  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar t(d),   g(d);   // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  t.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );

  t.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), +1.0));
  t.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), +1.0));

  t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Vector uvw(d);                          // velocity
  Enthalpy enth(t, g, uvw, time, solver, &fluid);
  enth.diffusion_set(TimeScheme::backward_euler());

  AC multigrid( &enth );

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    enth.new_time_step();
    // enth.smooth(0.0001, "temperature");
    multigrid.vcycle(ResRat(1e-4));

    if(time.current_step() % 10 == 0) {
      boil::plot->plot(t,   "t",  time.current_step());
    }

  }
  boil::plot->plot(t,   "t",  time.current_step());

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
