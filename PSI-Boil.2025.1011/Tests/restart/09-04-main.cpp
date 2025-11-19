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

  if(argc==1){
    boil::oout<<"One command line argument is required!"<<"\n";
    boil::oout<<"./Boil wmin (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

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
  time.print_time(false);


  Krylov * solver = new CG(d, Prec::di());

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr(p, f, uvw, time, solver, &fluid);
  Momentum ns( uvw, xyz, time, solver, &fluid);

  AC multigrid( &pr );

  Location loc("monitor", d, NX1, NY/2, NZ/2);

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  int ts=0;
  std::fstream input;
  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    uvw.  load("uvw",   ts);
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
    boil::oout<<"### RESTART ###\n";
    boil::oout<<"### time_step= "<<ts<<"\n";
  } else {
    boil::oout<<"### START FROM SCRATCH ###\n";
  }

  /*------------+
  |  time-loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    ns.cfl_max();

    ns.new_time_step();

    ns.solve(ResRat(1e-2));

    p = 0.0;

    multigrid.vcycle(ResRat(1e-2));

    ns.project(p);

    loc.print(uvw, Comp::v());
 
    if( time.current_step() % 100 == 0)
      boil::plot->plot(uvw,  p, "uvw-p",  time.current_step());

    /* output for restart */
    if( boil::timer.current_min() > (wmin-5)  // exit 5 min before wmin
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      std::fstream output;
      output << std::setprecision(16);
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }
  }

  boil::timer.stop();
  boil::timer.report();
}
