#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =   1.0;
const real LY =   0.25;

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

  Times time(20, 0.002); /* ndt, dt */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar h  (d), g  (d); // enthalpy    
  Scalar t  (d);         // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 0.0, 0.0, 1.0 ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  h.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  h.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  h.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 133780 ) ); // @ -30
  h.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), 284300 ) ); // @ ~30
  h.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  h.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  h = 0.5*(133780 + 284300);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  LookUpTable CO2("co2.txt");

  Matter fluid(d);
  fluid.mu (0.00010989);
  fluid.rho(961.94);
  fluid.cp (2254.3);

  fluid.variable(Set::rho()); /* fluid is variable and a function of t */
  fluid.variable(Set::mu());  /* fluid is variable and a function of t */

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr  (p,   f,   uvw, time, solver, &fluid);
  Momentum ns  (uvw, xyz,      time, solver, &fluid);
  Enthalpy enth(h,   g,   uvw, time, solver, &fluid);

  AC multigrid( &pr );

//uvw.load("uvw");

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    t.look_up(h, CO2, Column(2), Column(0));
    fluid.look_up(Set::rho(), t, CO2, Column(0), Column(3));
    fluid.look_up(Set::mu(),  t, CO2, Column(0), Column(6));

    enth.discretize();
    ns  .discretize();
    pr  .discretize();
    pr  .coarsen();

    enth.new_time_step();
    enth.solve(ResRat(0.0001), "temperature");

    ns.cfl_max();
    ns.new_time_step();

    ns.solve(ResRat(0.0001));

    for(int i=0; i<p.ni(); i++)
      for(int j=0; j<p.nj(); j++)
        for(int k=0; k<p.nk(); k++)
          p[i][j][k] = 0.0;
    
    multigrid.vcycle(ResRat(1e-3));
    p.exchange();
    ns.project(p);
    pr.update_rhs();
    uvw.exchange();

    if(time.current_step() % 100 == 0) {
      boil::plot->plot(uvw,   "uvw",   time.current_step());
      boil::plot->plot(p,h,t, "p-h-t", time.current_step());
      boil::plot->plot(*fluid.density(), "rho",  time.current_step());
    }

    OPR(p.min());
    OPR(p.max());
    OPR(h.min());
    OPR(h.max());
    OPR(t.min());
    OPR(t.max());
  }

  boil::plot->plot(uvw,   "uvw",   time.current_step()-1);
  boil::plot->plot(p,h,t, "p-h-t", time.current_step()-1);
  boil::plot->plot(*fluid.density(), "rho",  time.current_step()-1);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(p, "test", 0);
}
