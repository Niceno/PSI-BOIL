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
  |  physical properties  |  //four-fields-new-main.cpp
  +----------------------*/
  Matter fluid(d);
  fluid.rho   (1000.0);
  fluid.mu    (0.001); 

  Matter air(d);
  air.rho   (1.0);
  air.mu    (1.80e-5);

  Matter mixed(fluid, air, &c);
  mixed.sigma(0.072);

  Matter lagparticle(d);
  lagparticle.rho(1000.0);

  /*------------+
  |  time step  |
  +------------*/
  Times time(1000000000, 20.e-6);  /* ndt, dt */

  /*------------------+
  |  call Lagrangian  |
  +------------------*/
  Lagrangian lagp (c, &c, 0.56789, c, uvw, time, &mixed, &lagparticle);

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

  /*--------------------+
  |  initial condition  |  //main-bubble-3d.cpp--init-0
  +--------------------
  const Comp m = Comp::u();  //u0=1.0--z1
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k] = 1.0;*/
  //boil::plot->plot(uvw, p, c_sep, c_air, c_fluid, "initial", 0);
  for_vijk(c,i,j,k)
    c[i][j][k] = 0.0;

  /*----------+
  |  restart  |
  +----------*/
  uvw.load("uvw", 25600);
  t.  load("t",   25600);
  p.  load("p",   25600);
  //time.first_step(0);
  //boil::plot->plot(uvw, p, c_sep, c_air, c_fluid, "initial", 0);

  //head-format
  std::ofstream outfile("particles-trajectory.plt", std::ios_base::app);
  outfile << " VARIABLES= X, Y, Z, Up, Vp, Wp, Dp, ts " << boil::endl;
  outfile << " ZONE I= number  J= 1  K= 1" << boil::endl;
  outfile.close();

  /*---------------------+
  |  particle injection  |
  +---------------------*/
  int InjectionGroupNum = 100;
  for(int i=1; i<=InjectionGroupNum; i++) {
    boil::oout<< "num...i " << i << boil::endl;
    real xp_in_1 = 0.0005;
    real xp_in_2 = 0.0045;
    real xp_in_i = (xp_in_1-real((xp_in_2-xp_in_1)/real(InjectionGroupNum-1.0))) + real(i*(xp_in_2-xp_in_1)/real(InjectionGroupNum-1.0));
    //boil::oout<< "xp_in_i " << xp_in_i << boil::endl;

    lagp.add(Particle( Position(xp_in_i, 0.059999, 0.000), 
                       Diameter(10.e-6),
                       Position(0.0, -1.0, 0.0)));
  }

  /*#if resume    //bubble_colume_mini--main-aw.cpp
    //set resume to true when you want to resume the simulation
    //the numbers for csint, cbint, bck_ind: you find them in the file output
    disp.load("particles",12);
  #endif*/

  /*--------------------------------------+
  |  time loop part one : only particles  |
  +--------------------------------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# DT:        " << time.dt() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "########################" << boil::endl;

    boil::oout << "# particle time loop #" << boil::endl;

    lagp.advance(& xyz);

    if(time.current_step() % 20 == 0) {
      lagp.save("particles-trajectory", time.current_step());
    }
  }

  boil::oout << "particle-tracking-finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();
  getchar();

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
 '$Id: TJunction-TwoPhase.cpp,v 1.0 2018/11/22 13:13:00 MinZhu Exp $'/
+-----------------------------------------------------------------------------*/
