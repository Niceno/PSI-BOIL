#include "Include/psi-boil.h"

#define DIR 1

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
  #if DIR==0
    Body cyl("cylinder_x_128.stl");
  #elif DIR==1
    Body cyl("cylinder_y_128.stl");
  #else
    boil::oout << "Wrong direction!" << boil::endl;
    exit(0);
  #endif
  boil::plot->plot(cyl, "cylinder");

  /*-----------------------+
  |  computational domain  |
  +-----------------------*/
  Grid1D gx( Range<real>(-0.1001, 0.1001),   4, Periodic::yes() ); 
  Grid1D gy( Range<real>(-4.0001, 4.0001), 256, Periodic::yes() ); 
  Grid1D gz( Range<real>(-4.0004,12.0004), 512, Periodic::no() ); 

  /* small
  Grid1D gy( Range<real>(-2.0001, 2.0001), 128, Periodic::yes() );
  Grid1D gz( Range<real>(-2.0004, 6.0004), 256, Periodic::no() );
  */

  #if DIR==0
    Domain d(gx, gy, gz, & cyl);
  #else
    Domain d(gy, gx, gz, & cyl);
  #endif

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d), solid(d);
  fluid.mu(0.01);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Times time(20, 0.005); /* ndt, dt */
	
  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar p(d),   g(d);   // temperature
  Vector uvw(d), xyz(d); // velocity    

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet(), 0.,0.,1.) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  }

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr(p, g, uvw, time, solver, &fluid);
  Momentum ns(uvw, xyz, time, solver, &fluid);

  AC multigrid( &pr );

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    ns.cfl_max();
    ns.new_time_step();

    ns.solve(ResRat(0.01));

    p = 0.0;
    multigrid.vcycle(ResRat(1e-2));

    ns.project(p);
    pr.update_rhs();

    if(time.current_step() % 5 == 0) 
      boil::plot->plot(uvw, p, "uvw-p",  time.current_step());
    
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
