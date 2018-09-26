#include "Include/psi-boil.h"

#include <vector>

/* boundary conditions */
const real mu     = 6.0e-5;  /* => Re_tau = 1000 */
const real p_drop = 0.0036;  /* => Re_tau = 1000 */ 
const real t_drop = 0.5;     /* an arbitrary value */  
const real Pr_t   = 0.9;     /* turbulent prandtl number */

const real LX =  12.0;
const real LY =   2.0;
const real LZ =   4.0;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  const real b_des = 1.0;
  real       b_new = 0;

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx(Range<real>( 0.0, LX ), 64, Periodic::yes());
  Grid1D gy(Range<real>( -LY/2.0,  LY/2.0), 32, Periodic::no());
  Grid1D gz(Range<real>( 0.0, LZ ), 64, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  Times time(100, 0.04); /* ndt, dt */
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu    (mu);
  fluid.lambda(1.0e-3);
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p(d),   f(d);   // pressure
  Scalar q(d);
  Scalar wd(d), ws(d);   // wall distance and "source"
  Scalar mu_t(d);        // turbulent viscosity
  Scalar diff_t(d);      // turbulent diffusivity
  Scalar t(d),   g(d);   // temperature

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
  
  t.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::neumann(), t_drop ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::neumann(), t_drop ) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  wd.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  wd.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  wd.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), 0.0 ) );
  wd.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), 0.0 ) );
  wd.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  wd.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

  mu_t   = wd.shape();
  diff_t = wd.shape();

  Momentum ns(uvw, xyz,      time, solver, & fluid);
  Pressure pr(p,   f,   uvw, time, solver, & fluid);
  Enthalpy en(t,   g,   uvw, time, solver, & fluid);
  Distance di(wd,  ws,  uvw, time, solver);
  ns.convection_set(ConvScheme::central());

  AC multigrid( &pr );
  multigrid.stop_if_diverging( false );

  OMS(A);
  di.compute();
  OMS(B);
  boil::plot->plot(wd, "wd");

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

  /*--------+
  |  model  |
  +--------*/
  Model tm;

  /*------------+
  |  time loop  |
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

    const Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = p_drop * uvw.dV(m,i,j,k);

    for_vijk(t,i,j,k)
      g[i][j][k] = - t_drop * t.dV(i,j,k);

    // tm.smagorinsky( & ns, & mu_t, 0.065 );
    tm.wale( & ns, & mu_t, 0.1 );
 
    diff_t = mu_t / Pr_t;
    en.discretize( & diff_t ); 
    en.new_time_step();
    en.solve(ResRat(0.01), "enthalpy");

    ns.discretize( & mu_t ); 
    tm.tau_wall( & ns, wd, & xyz );

    ns.new_time_step();

    ns.solve(ResRat(0.01));

    p = 0.0;
    
    multigrid.vcycle(ResRat(0.01));

    p.exchange();   // is this needed? check!
    ns.project(p);
    pr.update_rhs();

    OPR(p.min());
    OPR(p.max());

    /*-------+
    |  plot  |
    +-------*/
    if( time.current_step() % 10000 == 0 ) {
      boil::plot->plot(mu_t,       "mu_t",   time.current_step());
      boil::plot->plot(uvw, p, t, "uvw-p-t", time.current_step());
    }

    /*-------+
    |  save  |
    +-------*/
    if( time.current_step() % 5000 == 0 ) {
      uvw.save ("uvw",  time.current_step());
      mu_t.save("mu_t", time.current_step());
      p.save   ("p",    time.current_step());
    }
  }

  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: main-turbulent-heat.cpp,v 1.1 2011/05/27 21:30:59 niceno Exp $'/
+-----------------------------------------------------------------------------*/
