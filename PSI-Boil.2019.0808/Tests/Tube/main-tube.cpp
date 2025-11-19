#include "Include/psi-boil.h"

#define DIR 0

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*-------------------------+
  |  specify flow direction  |
  +-------------------------*/
  #if DIR==0 
    Body cyl("tube_x.stl");
  #endif
  #if DIR==1 
    Body cyl("tube_y.stl");
  #endif
  #if DIR==2 
    Body cyl("tube_z.stl");
  #endif

  /*----------------+
  |  immersed body  |
  +----------------*/
  boil::plot->plot(cyl, "tube");

  /*-----------------------+
  |  computational domain  |
  +-----------------------*/
  Grid1D gx( Range<real>(-0.75001, 0.75001),  32, Periodic::yes() ); 

  Domain d(gx, gx, gx, & cyl);

  /*----------------+
  |  immersed body  |
  +----------------*/
  boil::plot->plot(cyl, "tube-2");

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d), solid(d);
  fluid.mu    (0.01);
  fluid.lambda(0.01);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Times time(100, 0.005); /* ndt, dt */
	
  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar p(d),   g(d);   // temperature
  Vector uvw(d), xyz(d); // velocity    

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr(p, g, uvw, time, solver, &fluid);
  Momentum ns(uvw, xyz, time, solver, &fluid);

  AC multigrid( &pr );

  const real b_des = 1.0;
  real       b_new = 0;

  Comp m;
  if(DIR==0) m=Comp::u();
  if(DIR==1) m=Comp::v();
  if(DIR==2) m=Comp::w();
  for_vmijk(xyz,m,i,j,k)
    xyz[m][i][j][k] = 1.0 * uvw.dV(m,i,j,k);

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
/*
    p = 0.0;
    multigrid.vcycle(ResRat(1e-2));

    ns.project(p);
    pr.update_rhs();

    if(time.current_step() % 10 == 0) 
      boil::plot->plot(uvw, p, "uvw-p",  time.current_step());
*/    
    /*--------------------------+ 
    |  recompute pressure drop  |
    +--------------------------*/ 
    b_new = ns.bulk(m, 0.0);
    real p_drop = fluid.rho() * (b_des - b_new) / time.dt();

    OPR(b_des);
    OPR(b_new);

    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = p_drop * uvw.dV(m,i,j,k);

  }

  boil::plot->plot(uvw, p, "uvw-p",  time.current_step()-1);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}
