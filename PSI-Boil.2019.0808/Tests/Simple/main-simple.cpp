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

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu(0.01);
  fluid.rho(2.0);

  Times time(10, 0.002); /* ndt, dt */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar press(d);

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
  
  press.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &fluid);
  ns.convection_set(TimeScheme::backward_euler());
  ns.diffusion_set(TimeScheme::backward_euler());

  Pressure pr(p, f, uvw, time, solver, &fluid);
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
	  
    /* momentum */
    ns.cfl_max();
    ns.new_time_step();

    // NEW: SIMPLE ITERATIONS
    for(int in=0; in<6; in++) {
      ns.grad(press);
      ns.convection();
      ns.solve(ResRat(0.0001));

      for_vijk(p,i,j,k) p[i][j][k] = 0.0;
      p.exchange();

      multigrid.vcycle(ResRat(1e-2));
      p.exchange();
      ns.project(p);

      for_vijk(p,i,j,k)
        press[i][j][k] += p[i][j][k];
      press.exchange();
    }

    if(time.current_step() % 100 == 0) {
      boil::plot->plot(uvw,"uvw", time.current_step());
      boil::plot->plot(p,   "p",  time.current_step());
    }

  }
  boil::plot->plot(uvw,  "uvw",  time.current_step());
  boil::plot->plot(p,    "p",    time.current_step());
  boil::plot->plot(press,"press",time.current_step());

  //uvw.plot_par("uvw", 10);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(p, "test", 0);
}
