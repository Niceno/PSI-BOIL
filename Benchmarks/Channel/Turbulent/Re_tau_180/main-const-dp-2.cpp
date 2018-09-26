#include "Include/psi-boil.h"

#include <vector>

/* boundary conditions */
const real mu =   0.0003333333333333333333333333333;

/******************************************************************************/
main(int argc, char * argv[]) {

  const real b_des = 1.0;
  real       b_new = 0;

  #include "domain-2.cpp"

  Times time(40000, 0.015); /* ndt, dt */
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu(mu);
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p(d),   f(d);   // pressure
  Scalar q(d);

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  Momentum ns(uvw, xyz,      time, solver, & fluid);
  Pressure pr(p,   f,   uvw, time, solver, & fluid);
  ns.convection_set(ConvScheme::central());

  AC multigrid( &pr );

  /*--------------+
  |  random flow  |
  +--------------*/
  int ne = 100;
  real tt = 0.1; // turbulent time scale
  real tl = 0.1; // turbulent length scale
  RandomFlow rf( ne, tt, tl );

  uvw.randomize(rf, 0.5, 0.0);

  const Comp m = Comp::u();
  for_vmijk(xyz,m,i,j,k)
    uvw[m][i][j][k] += 1.0;

  boil::plot->plot(uvw,"uvw", 0);

  /*------------+
  |  time loop  |
  +------------*/
  boil::timer.start();
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

    ns.solve(ResRat(0.0001));

    p = 0.0;
    
    if( multigrid.vcycle(ResRat(0.01)) ) {
      p.exchange();   // is this needed? check!
      ns.project(p);
      uvw.exchange(); // this is not needed. check
      pr.update_rhs();
    }

    OPR(p.min());
    OPR(p.max());

    /*---------+
    |  p_drop  |
    +---------*/
    b_new = ns.bulk(Comp::u(), LX*0.333);
    real p_drop = 0.0036;

    const Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = p_drop * uvw.dV(m,i,j,k);

    OPR( b_new );
    OPR( p_drop );

    /*-------+
    |  plot  |
    +-------*/
    if( time.current_step() % 4000 == 0 ) {
      ns.get_q(&q);
      boil::plot->plot(q,  "q",   time.current_step());
      boil::plot->plot(p,  "p",   time.current_step());
      boil::plot->plot(uvw,"uvw", time.current_step());
    }
    /*-------+
    |  save  |
    +-------*/
    if( time.current_step() % 100 == 0 ) {
      uvw.save("uvw", time.current_step());
    }
  }
  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: main-const-dp-2.cpp,v 1.6 2011/05/25 11:09:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
