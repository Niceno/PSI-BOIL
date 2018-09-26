#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =   1.0;
const real LY =   0.2;
const int  NX =  64;
const int  NY =  32;
const int  NZ =  32;

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
  Grid1D g_x( Range<real>(0, LX), 
              NX, Periodic::no());

  Grid1D g_y( Range<real>(-0.5*LY, 0.5*LY),  
              NY, Periodic::yes());

  Grid1D g_z( Range<real>(0,LY),  
              Range<real>( LY/128.0, LY/128.0 ),
              NZ, Periodic::no());

  /*----------+
  |  domains  |
  +----------*/
  Domain dom_1(g_x, g_y, g_z, "dom_long_x");
  Domain dom_2(g_x, g_y, g_z, "dom_long_y");

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid_1(dom_1);
  Matter fluid_2(dom_2);
  fluid_1.mu    (0.001);
  fluid_2.mu    (0.001);
  fluid_1.lambda(0.001);
  fluid_2.lambda(0.001);

  Times time(500, 0.0025); /* ndt, dt */
	
  /*-----------------+
  |  linear solvers  |
  +-----------------*/
  Krylov * solver_1 = new CG(dom_1);
  Krylov * solver_2 = new CG(dom_2);

  /*-----------------------------------+
  |  define unknowns and default b.c.  |
  +-----------------------------------*/
  Vector uvw_1(dom_1), xyz_1(dom_1); // vel 1
  Vector uvw_2(dom_2), xyz_2(dom_2); // vel 2
  Scalar p_1  (dom_1), f_1  (dom_1); // p 1
  Scalar p_2  (dom_2), f_2  (dom_2); // p 2
  Scalar t_1  (dom_1, "temperature_1"), g_1  (dom_1); // t 1
  Scalar t_2  (dom_2, "temperature_2"), g_2  (dom_2); // t 2

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw_1.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 1.0, 0.0, 0.0 ) );
    uvw_1.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw_1.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw_1.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw_1.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw_1.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );

    uvw_2.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 1.0, 0.0, 0.0 ) );
    uvw_2.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw_2.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw_2.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw_2.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw_2.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }  

  p_1.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p_1.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p_1.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p_1.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  p_2.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p_2.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p_2.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p_2.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  t_1.bc().add( BndCnd( Dir::imin(), BndType::inlet(), "0.25*(1.0 - (y*10.0)^2)" ) );
  t_1.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  t_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  t_1.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 1.0 ) );
  t_1.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );

  t_2.bc().add( BndCnd( Dir::imin(), BndType::inlet(), 0.0 ) );
  t_2.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  t_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  t_2.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 1.0 ) );
  t_2.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );

  real ** copy_plane;

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr_1( p_1,   f_1,   uvw_1, time, solver_1, &fluid_1 );
  Pressure pr_2( p_2,   f_2,   uvw_2, time, solver_2, &fluid_2 );
  Momentum ns_1( uvw_1, xyz_1,        time, solver_1, &fluid_1 );
  Momentum ns_2( uvw_2, xyz_2,        time, solver_2, &fluid_2 );
  Enthalpy en_1( t_1,   g_1,   uvw_1, time, solver_1, &fluid_1 );
  Enthalpy en_2( t_2,   g_2,   uvw_2, time, solver_2, &fluid_1 );

  AC sol_1( &pr_1 );
  AC sol_2( &pr_2 );

  sol_1.stop_if_diverging(false);
  sol_2.stop_if_diverging(false);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "domain 1" << boil::endl;

    en_1.discretize();
    en_1.new_time_step();
    en_1.solve(ResRat(0.01), "temperature 1");
    
    ns_1.cfl_max();
    ns_1.new_time_step();

    ns_1.solve(ResRat(0.01));

    for(int i=0; i<p_1.ni(); i++)
      for(int j=0; j<p_1.nj(); j++)
        for(int k=0; k<p_1.nk(); k++)
          p_1[i][j][k] = 0.0;
    
    sol_1.vcycle(ResRat(0.01));
    p_1.exchange();
    ns_1.project(p_1);
    uvw_1.exchange();

    boil::oout << "domain 2" << boil::endl;

    /* copy and paste the boundary */
    t_1.bnd_extract( Dir::imax(), &copy_plane );
    t_2.bnd_insert ( Dir::imin(),  copy_plane );

    en_2.discretize();
    en_2.new_time_step();
    en_2.solve(ResRat(0.01), "temperature 2");

    ns_2.cfl_max();
    ns_2.new_time_step();

    ns_2.solve(ResRat(0.01));

    for(int i=0; i<p_2.ni(); i++)
      for(int j=0; j<p_2.nj(); j++)
        for(int k=0; k<p_2.nk(); k++)
          p_2[i][j][k] = 0.0;
    
    sol_2.vcycle(ResRat(0.01));
    p_2.exchange();
    ns_2.project(p_2);
    uvw_2.exchange();

    if(time.current_step() % 100 == 0) {
      boil::plot->plot(uvw_1,"uvw_1", time.current_step());
      boil::plot->plot(uvw_2,"uvw_2", time.current_step());
      boil::plot->plot(t_1,  "t_1",   time.current_step());
      boil::plot->plot(t_2,  "t_2",   time.current_step());
    }
  }

  //uvw.plot_par("uvw", 10);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw_2, p_2, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-two-domains.cpp,v 1.20 2014/10/15 13:44:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
