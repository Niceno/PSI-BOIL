#include "Include/psi-boil.h"

#include <vector>

/* parameters */
const real LX  =  0.06;
const real LY  =  0.051;
const real rho =  1.205;
const real nu  = 15.11e-4; /* realistic is 15.11e-6 */
const real mu  = rho * nu; 

/****************************************************************************/
main(int argc, char * argv[]) {

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  boil::timer.start();

  const real b_des = 3.86;
  real       b_new = 0;

  Times time(6000, 0.00004); /* ndt, dt */

  /*--------+
  |  grids  |
  +--------*/
  Grid1D gx(Range<real>( -LX/2.0,  LX/2.0 ), 80, Periodic::yes());
  Grid1D gy(Range<real>( 0, LY ), 64, Periodic::no());

  /*----------------+
  |  immersed body  |
  +----------------*/
  Body cube("cube.stl");
  boil::plot->plot(cube, "cube");

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gy, & cube);
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d), solid(d);
  fluid.mu (mu);
  fluid.rho(rho);
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p  (d), f  (d); // pressure
  Scalar q  (d);         // q

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall(), 0.0, 0.0, 0.0 ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall(), 0.0, 0.0, 0.0 ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  
  const Comp m = Comp::u();
  for_vmijk(xyz,m,i,j,k)
    xyz[m][i][j][k] = 1.0 * uvw.dV(m,i,j,k);

  Pressure pr( p,   f,   uvw, time, solver, &fluid );
  Momentum ns( uvw, xyz,      time, solver, &fluid );
  ns.convection_set(ConvScheme::central());

  AC multigrid( &pr );

  for(time.start(); time.end(); time.increase()) {

    ns.discretize();
    ns.cfl_max();
    ns.new_time_step(); 

    ns.solve(ResRat(1e-2));

    p = 0.0;
    
    multigrid.vcycle(ResRat(0.01));

    ns.project(p);

    pr.update_rhs();

    b_new = ns.bulk(Comp::u(), LX*0.333);

    real p_drop = fluid.rho()->value() * (b_des - b_new) / time.dt();

    const Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = p_drop * uvw.dV(m,i,j,k);

    if(time.current_step() % 120 == 0) {
      boil::plot->plot(uvw, p, "uvw,p",  time.current_step());
    }

    OPR( b_new );
    OPR( p_drop );
  }
  boil::plot->plot(uvw,p,"uvw,p", time.current_step()-1);

  boil::timer.stop();
  boil::timer.report();
}
