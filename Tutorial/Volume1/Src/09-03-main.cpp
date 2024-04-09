#include "Include/psi-boil.h"

#include <vector>

/****************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*---------------------------------+
  |  1. grids, obstacles and domain  |
  +---------------------------------*/
  #include "09-03-common.h"
	
  /* timer */
  Times time(100, 0.00002); /* ndt, dt */

  /*---------------------+
  |  2. define unknowns  |
  +---------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p  (d), f  (d); // pressure

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
  
  /*------------------------------------+
  |  3. physical properties and solver  |
  +------------------------------------*/
  Matter fluid(d);
  fluid.mu (mu);
  fluid.rho(rho);
	
  Krylov * solver = new CG(d, Prec::di());

  /*-------------------------+
  |  4. transport equations  |
  +-------------------------*/
  Pressure pr( p,   f,   uvw, time, solver, &fluid );
  Momentum ns( uvw, xyz,      time, solver, &fluid );
  ns.convection_set(ConvScheme::central());

  /*-----------------------------------+
  |  5. preparation for the time loop  |
  +-----------------------------------*/
  AC multigrid( &pr );

  Comp m = Comp::u();
  for_vmijk(xyz,m,i,j,k) xyz[m][i][j][k] = 1.0 * uvw.dV(m,i,j,k);

  /*---------------+
  |  6. time loop  |
  +---------------*/
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

    ns.solve(ResRat(0.001));

    p = 0.0;
    
    multigrid.vcycle(ResRat(0.001));

    ns.project(p);

    pr.update_rhs();

    real b_new = ns.bulk(Comp::u(), LX*0.333);

    real p_drop = fluid.rho() * (b_des - b_new) / time.dt();

    Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k) xyz[m][i][j][k] = p_drop * uvw.dV(m,i,j,k);

    OPR( b_new );
    OPR( p_drop );

    /* save */
    if( time.current_step() % 5 == 0 ) uvw.save("uvw", time.current_step());
    if( time.current_step() % 5 == 0) 
      boil::plot->plot(uvw,p,"uvw,p", time.current_step());
  }

  boil::timer.stop();
  boil::timer.report();
}
