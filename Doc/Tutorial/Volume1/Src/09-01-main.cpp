#include "Include/psi-boil.h"

/* dimensions */
const real H  = 0.41;
const real L  = 2.2;

/* resolutions */
const int NY =  64;
const int NX1 = NY;
const int NX2 = NY;
const int NZ  = NY/8;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  boil::plot = new PlotTEC();

  /*----------------+
  |  immersed body  |
  +----------------*/
  Body cyl("09-01-cylinder.stl");

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gz ( Range<real>(-H/4, H/4), NZ,  Periodic::yes());
  Grid1D gy ( Range<real>( 0.0, H  ), NY,  Periodic::no());
  Grid1D gx1( Range<real>( 0.0, H  ), NX1, Periodic::no());
  const real dx = H/(real)NY;
  Grid1D gx2( Range<real>( H,   L  ), Range<real>(dx, 8.0*dx), NX2, Periodic::no());
  Grid1D gx( gx1, gx2, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz, &cyl);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); /* velocity and its source */ 
  Scalar p  (d), f  (d); /* pressure and its source */ 

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 
                   "4.0*1.5*y*(0.41-y)/0.41^2", "0.0", "0.0") );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);

  fluid.mu(0.001);

  Times time(1000, dx/8.0); 

  Krylov * solver = new CG(d, Prec::di());

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr(p, f, uvw, time, solver, &fluid);
  Momentum ns( uvw, xyz, time, solver, &fluid);

  AC multigrid( &pr );

  Location loc("monitor", d, NX1, NY/2, NZ/2);

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
	  
    ns.cfl_max();

    ns.new_time_step();

    ns.solve(ResRat(1e-2));

    p = 0.0;

    multigrid.vcycle(ResRat(1e-2));

    ns.project(p);

    loc.print(uvw, Comp::v());
 
    if( time.current_step() % 100 == 0)
      boil::plot->plot(uvw,  p, "uvw,p",  time.current_step());
  }

  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: 09-01-main.cpp,v 1.7 2011/05/25 11:10:29 niceno Exp $'/
+-----------------------------------------------------------------------------*/
