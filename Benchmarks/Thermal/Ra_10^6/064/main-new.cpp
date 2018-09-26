#include "Include/psi-boil.h"

/* parameters */
const real LX =   1.0;
const real LY =   0.125;

const int NX = 64;
const int NY =  4;

const real Pr = 0.71;
const real Ra = 1.0e+6;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

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
  fluid.mu( Pr );

  Times time(12000, 0.000025); /* ndt, dt */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar t  (d), g  (d); // t.

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
  
  t.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), +0.5 ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), -0.5 ) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  Location m0("m0", d, NX/2, NY/2, NX/4);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr  ( p,   f,   uvw, time, solver, &fluid);
  Momentum ns  ( uvw, xyz,      time, solver, &fluid);
  Enthalpy enth( t,   g,   uvw, time, solver, &fluid);

//  ns.diffusion (TimeScheme::backward_euler());
//  ns.convection(TimeScheme::backward_euler());

  AC multigrid( &pr );

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    enth.new_time_step();
    enth.convection();
    enth.solve(ResRat(0.0001));

    ns.cfl_max();
    ns.new_time_step();

    for_vij(g,i,j)
      for(int k=g.si()+1; k<=g.ei(); k++)
         xyz[Comp::w()][i][j][k] = 0.5 * Pr * Ra * (t[i][j][k] + t[i][j][k-1]);

    ns.convection();
    ns.solve(ResRat(0.0001));

    for_avijk(p,i,j,k) p[i][j][k] = 0.0;

    multigrid.vcycle(ResRat(0.001));
    p.exchange();
    ns.project(p);

    pr.update_rhs();

    m0.print(uvw,Comp::u());
  }

  boil::plot = new PlotGMV();
  boil::plot->plot(uvw, t, "uvw,t", time.current_step()-1);
  boil::plot->plot(p,      "p",     time.current_step()-1);

  Rack    r0("u-comp", d, NX/2+1, NY/2, Range<int>(1,NX));
  Rack    r2("w-comp", d, Range<int>(1,NX), NY/2, NX/2+1);
  r0.print(uvw,Comp::u());
  r2.print(uvw,Comp::w());

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: main-new.cpp,v 1.11 2011/05/25 11:09:47 niceno Exp $'/
+-----------------------------------------------------------------------------*/
