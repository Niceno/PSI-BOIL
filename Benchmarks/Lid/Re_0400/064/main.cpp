#include "Include/psi-boil.h"

/* parameters */
const real LX =   1.0;
const real LY =   0.0125;

const int NX = 64;
const int NY =  4;

const real Re = 400.0;

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
  Grid1D gx( Range<real>( 0,LX ), 
             Range<real>( LX/NX, LX/NX ),
              NX, Periodic::no());

  Grid1D gy( Range<real>( 0,LY ),  
             Range<real>( LY/NY, LY/NY ),
              NY, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu(1.0/Re);

  Times time(3200, 0.025); /* ndt, dt */
	
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
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::inlet(), 1.0, 0.0, 0.0 ) );
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
  
  Location m0("m0", d, NX/2, NY/2, NX/4);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr( p,   f,   uvw, time, solver, &fluid);
  Momentum ns( uvw, xyz,      time, solver, &fluid);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(TimeScheme::backward_euler());

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

    // NEW: SIMPLE ITERATIONS
    for(int in=0; in<6; in++) {
      ns.grad(press);
      ns.convection();
      ns.solve(ResRat(0.0001));

      for_avijk(p,i,j,k) p[i][j][k] = 0.0;

      multigrid.vcycle(ResRat(0.001));
      p.exchange();
      ns.project(p);

      for_vijk(p,i,j,k)
        press[i][j][k] += p[i][j][k];
      press.exchange();

      uvw.exchange(); // to be sure
      if( pr.update_rhs() < 1.0e-8) break;
    }

    m0.print(uvw,Comp::u());
  }
  boil::plot->plot(uvw,  "uvw",  time.current_step()-1);
  boil::plot->plot(p,    "p",    time.current_step()-1);
  boil::plot->plot(press,"press",time.current_step()-1);

  Rack    r0("u-comp", d, NX/2+1, NY/2, Range<int>(1,NX));
  Rack    r2("w-comp", d, Range<int>(1,NX), NY/2, NX/2+1);
  r0.print(uvw,Comp::u());
  r2.print(uvw,Comp::w());

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.11 2011/05/25 11:09:37 niceno Exp $'/
+-----------------------------------------------------------------------------*/
