#include "Include/psi-boil.h"

/* parameters */
const real LX =   1.0;
const real LY =   0.125;

const int NX = 64;
const int NY =  4;

const real Pr = 0.71;
const real Ra = 1.0e+5;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>( -0.5*LX, 0.5*LX ), NX, Periodic::no());
  Grid1D gy( Range<real>( 0, LY ),           NY, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);

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
  
  /*---------------------------------------+
  |  physical properties, time and solver  |
  +---------------------------------------*/
  Matter fluid(d);
  fluid.mu( Pr );

  Times time(4000, 0.000075); /* ndt, dt */
	
  Krylov * solver = new CG(d, Prec::di());

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr  ( p,   f,   uvw, time, solver, &fluid);
  Momentum ns  ( uvw, xyz,      time, solver, &fluid);
  Enthalpy enth( t,   g,   uvw, time, solver, &fluid);

  AC multigrid( &pr );

  Location m0("m0", d, t.sI()+NX/2, t.sJ()+NY/2, t.sK()+NX/4);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    enth.new_time_step();
    enth.solve(ResRat(0.001));

    ns.cfl_max();
    ns.new_time_step();

    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = Pr*Ra * 0.5*(t[i][j][k]+t[i][j][k-1]) * xyz.dV(m,i,j,k);

    ns.solve(ResRat(0.001));

    multigrid.vcycle(ResRat(0.001));
    p.exchange();
    ns.project(p);

    pr.update_rhs();

    m0.print(uvw,Comp::u());
  }

  boil::plot = new PlotTEC();
  boil::plot->plot(uvw, t, "uvw-t", time.current_step()-1);
  boil::plot->plot(p,      "p",     time.current_step()-1);

#if 0
  Rack    r0("u-comp", d, t.sI()+NX/2+1, t.sJ()+NY/2, Range<int>(t.sK(),t.sK()+NX));
  Rack    r2("w-comp", d, Range<int>(t.sI(),t.sI()+NX), t.sJ()+NY/2, t.sK()+NX/2+1);
  r0.print(uvw,Comp::u());
  r2.print(uvw,Comp::w());
#endif

  const int ii = t.si();   /* i inside the domain */
  const int iw = ii-1;     /* i in the wall */
  const int j  = t.ej()/2; /* j in the middle */
  boil::oout << " Nusselt number " << boil::endl;
  for_vk(t,k) {
    real nu = (t[iw][j][k] - t[ii][j][k]) / t.dxw(ii);
    boil::oout << t.zc(k) << " " << nu << boil::endl;
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}
