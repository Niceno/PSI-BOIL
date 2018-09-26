#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =   1.2;
const real LY =   0.05;
const real LZ =   0.14;

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
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 
             Range<real>( LX/256.0, LX/256.0 ),
             256, Periodic::no());

  Grid1D gy( Range<real>(0,LY),  
             Range<real>( LY/ 4.0, LY/ 4.0 ),
              4, Periodic::yes());

  Grid1D gz( Range<real>(0,LZ), 
             Range<real>( LZ/ 32.0, LZ/ 32.0 ),
              32, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  // orig Times time(320, 0.002); /* ndt, dt */
  Times time(400, 0.001); /* ndt, dt */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d); // concentration
  Scalar press(d);

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
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
  
  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  for_vijk(c,i,j,k) 
    c[i][j][k] = 1.0;

  for_vijk(c,i,j,k) {
    if(c.xc(i) < 0.0 && c.zc(k) < 0.1 )
      c[i][j][k] = 0.0;
  }
  
  c.exchange();
  boil::plot->plot(c, "c0", 0);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (0.00001 * 1.0e+1);
  air  .rho   (1.);
  water.mu    (0.001 * 1.0e+1);
  water.rho   (1000.0);

  Matter mixed(air, water, c);

  ColorFunction conc  (c,   g, 1.0, 1.0, uvw, time, solver); 
  conc.sharpen();
  boil::plot->plot(c, "c1", 0);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::upwind());

  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging( false );

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    conc.new_time_step();
    OMS("calling_convection");
    conc.convection();
    conc.advance();
    conc.sharpen();

    const Comp m = Comp::w();
    for(int i=1; i<xyz.ni(m)-1; i++)
      for(int j=1; j<xyz.nj(m)-1; j++)
        for(int k=1; k<xyz.nk(m)-1; k++)
          xyz[m][i][j][k] = -9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();
    ns.convection();
    ns.solve(ResRat(0.0001));

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(0.01));
    p.exchange();
    ns.project(p);

    pr.update_rhs();
    ns.cfl_max();

    if(time.current_step() % 100 == 0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
      boil::plot->plot(p,      "p",      time.current_step());
      boil::plot->plot(press,  "press",  time.current_step());
    }

  }

  //uvw.plot_par("uvw", 10);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw,c, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-dam.cpp,v 1.1 2015/02/19 09:28:34 sato Exp $'/
+-----------------------------------------------------------------------------*/
