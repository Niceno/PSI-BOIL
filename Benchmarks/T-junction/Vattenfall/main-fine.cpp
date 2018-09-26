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
  Body pipe("vattenfall_B.stl");

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx ( Range<real>(-0.2399,  1.0399), 512, Periodic::no());
  Grid1D gy ( Range<real>(-0.2399,  0.0799), 128, Periodic::no());
  Grid1D gz ( Range<real>(-0.0799,  0.0799),  64, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz, &pipe);

  boil::plot->plot(pipe, "pipe");

  Times time(8000, 0.0005); /* ndt, dt */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);  // velocity
  Scalar p  (d), f  (d);  // pressure
  Scalar t  (d), g  (d);  // temperature
  Scalar mu_t(d);

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 0.39, 0.00, 0.00) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::inlet(), 0.00, 0.38, 0.00) );
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
  
  t.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 15.0) );
  t.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), 30.0) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu    (   0.001);
  fluid.lambda(   0.01 );
  fluid.rho   (1000.0  );
  fluid.cp    (4183.0  );

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &fluid);
  ns.convection_set(ConvScheme::central());
  Enthalpy enth(t,   g,   uvw, time, solver, &fluid);
  Pressure pr  (p,   f,   uvw, time, solver, &fluid);

  AC multigrid( &pr );

//  Location loc_1("monitor_1", d, d.I(0.35), NY/2, NZ/2);
//  Location loc_2("monitor_2", d, d.I(0.45), NY/2, NZ/2);
//  Location loc_3("monitor_3", d, d.I(0.55), NY/2, NZ/2);
//  Location loc_4("monitor_4", d, d.I(0.65), NY/2, NZ/2);
//  Location loc_5("monitor_5", d, d.I(0.75), NY/2, NZ/2);
 
  /*------------------+ 
  |  physical models  |
  +------------------*/
  Model mod(d);

  /*--------------------------+ 
  |  force profiles to inlet  |
  +--------------------------*/
  Profile u_in("vattenfall_in_U.pro");
  Profile v_in("vattenfall_in_U.pro");

  u_in.scale_coords(0.07); /* D = 0.14 */
  v_in.scale_coords(0.05); /* D = 0.1  */
  u_in.scale_values(0.5);  
  v_in.scale_values(0.5);  

  Comp m=Comp::u();
  for_vmjk(uvw, m, j, k) {
    real r = sqrt( uvw.yc(m,j)*uvw.yc(m,j) + uvw.zc(m,k)*uvw.zc(m,k) );
    uvw[m][uvw.si(m)-1][j][k] = u_in.value_at(r);
  }

  m=Comp::v();
  for_vmik(uvw, m, i, k) {
    real r = sqrt( uvw.xc(m,i)*uvw.xc(m,i) + uvw.zc(m,k)*uvw.zc(m,k) );
    uvw[m][i][uvw.sj(m)-1][k] = v_in.value_at(r);
  }

  /*----------+
  |  restart  |
  +----------*/
  uvw. load("uvw", 7500);
  t.   load("t",   7500);
  time.first_step (7500);

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
	  
    /*--------------------+
    |  smagorinsky model  |
    +--------------------*/
    mod.smagorinsky( &mu_t, uvw, fluid );

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

//    /*-----------------------------+
//    |  print monitoring locations  |
//    +-----------------------------*/
//    loc_1.print(uvw, 1);
//    loc_2.print(uvw, 1);
//    loc_3.print(uvw, 1);
//    loc_4.print(uvw, 1);
//    loc_5.print(uvw, 1);
//
    /*------------+
    |  save/plot  |
    +------------*/
    if(time.current_step() % 2000 == 0) {
      uvw.  save("uvw",  time.current_step());
      p.    save("p",    time.current_step());
      t.    save("t",    time.current_step());
    }

    if(time.current_step() % 500 == 0) {
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
 '$Id: main-fine.cpp,v 1.3 2011/05/25 11:09:47 niceno Exp $'/
+-----------------------------------------------------------------------------*/
