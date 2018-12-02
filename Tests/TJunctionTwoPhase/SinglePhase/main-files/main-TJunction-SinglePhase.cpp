#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------------+
  |  immersed body  |
  +----------------*/
  Body pipe("mins-triangles.stl");

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>( -0.1050, 0.1050), 21*64, Periodic::no());
  Grid1D gy( Range<real>(  0.0000, 0.0600),  6*64, Periodic::no());
  Grid1D gz( Range<real>( -0.0100, 0.0100),    2, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz, &pipe);

  boil::plot->plot(pipe, "pipe");

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);  // velocity
  Scalar p  (d), f  (d);  // pressure
  Scalar t  (d), g  (d);  // temperature
  Scalar mu_t       (d);
  Scalar c          (d);


  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::inlet(), 0.00,-1.00, 0.00) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  t.bc().add( BndCnd( Dir::imin(), BndType::outlet() ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), 30.0) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  c = p.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  /*Matter fluid(d);
  fluid.rho   (1000.0);
  fluid.mu    (0.001); */

  Matter air(d);
  air.rho   (1.0);
  air.mu    (1.80e-5);

  /*Matter mixed(fluid, air, &c);
  mixed.sigma(0.072);*/

  /*------------+
  |  time step  |
  +------------*/
  Times time(1000000000, 20.e-6);  /* ndt, dt */

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Momentum ns( uvw, xyz, time, solver, &air);
  //ns.convection_set(ConvScheme::central());
  ns.convection_set(ConvScheme::muscl());
  Enthalpy enth(t,   g,   uvw, time, solver, &air);
  Pressure pr  (p,   f,   uvw, time, solver, &air);

  AC multigrid( &pr );

  //multigrid.stop_if_diverging(false);
  //multigrid.max_cycles(15);

  /*----------+
  |  restart  |
  +----------*/
  //uvw.load("uvw", 25600);
  //t.  load("t",   25600);
  //p.  load("p",   25600);
  //time.first_step(25600);
  //boil::plot->plot(uvw, p, c_sep, c_air, c_fluid, "initial", 0);

  /*------------+
  |  time-loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {
    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    boil::oout << "# fluid time loop #" << boil::endl;

    ns.discretize( &mu_t );

    enth.new_time_step();
    enth.solve(ResRat(1e-2), "t");

    for_vijk(t,i,j,k) {
      if( d.ibody().on(i,j,k) ) {
        if( t[i][j][k] > 30.0 ) t[i][j][k] = 30.0;
        if( t[i][j][k] < 15.0 ) t[i][j][k] = 15.0;
      } 
    }

    ns.cfl_max();

    ns.new_time_step();
    ns.solve(ResRat(1e-2));

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-2));

    p.exchange();
    ns.project(p);
    pr.update_rhs();

    /*-------------+
    |  dt control  |
    +-------------*/
    real cfl_limit = 0.30;
    real cflmax = ns.cfl_max();
    time.control_dt(cflmax, cfl_limit, 0.001);

    /*------------+
    |  save/plot  |
    +------------*/
    if(time.current_step() % 1000 == 0) {
      uvw.  save("uvw",  time.current_step());
      p.    save("p",    time.current_step());
      t.    save("t",    time.current_step());
    }

    if(time.current_step() % 100 == 0) {
      boil::plot->plot(uvw,  p, "uvw,p",  time.current_step());
      boil::plot->plot(uvw,  t, "uvw,t",  time.current_step());
      boil::plot->plot(mu_t,     "mu_t",  time.current_step());
      boil::plot->plot(c,           "c",  time.current_step());
    }

    /*---------+
    |  exit ?  |
    +---------*/
    std::ifstream infile;
    infile.open ("stop.now", std::ifstream::in);
    if(infile.good()) exit(0);
    infile.close();
  }
  boil::plot->plot(uvw,  p, "uvw,p",  time.current_step()-1);
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
/*-----------------------------------------------------------------------------+
 '$Id: TJunction-SinglePhase.cpp,v 1.0 2018/11/22 13:13:00 MinZhu Exp $'/
+-----------------------------------------------------------------------------*/
