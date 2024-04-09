#include "Include/psi-boil.h"

/* boundary conditions */
const real L   = 1.0;
const real mu  = 0.01;
const real lam = 0.01;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D g_per( Range<real>(0,L),  5, Periodic::yes());
  Grid1D g_non( Range<real>(0,L), 11, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d_px(g_per, g_non, g_non);
  Domain d_py(g_non, g_per, g_non);
  Domain d_pz(g_non, g_non, g_per);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter f_x(d_px), f_y(d_py), f_z(d_pz);

  Times time(1, 1e+3); /* ndt, dt */
	
  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar t_x(d_px), g_x(d_px);   // temperature
  Scalar t_y(d_py), g_y(d_py);   // temperature
  Scalar t_z(d_pz), g_z(d_pz);   // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  t_x.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  t_x.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  t_y.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t_y.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  t_z.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  t_z.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

  /* x */
  t_x.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  t_x.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  t_x.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), -1.0));
  t_x.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), +1.0));
  
  /* y */
  t_y.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  t_y.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  t_y.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), -1.0));
  t_y.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), +1.0));

  /* z */
  t_z.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), -1.0 ) );
  t_z.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), -1.0) );
  t_z.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), +1.0));
  t_z.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), +1.0));

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Vector u_x(d_px), 
         u_y(d_py), 
         u_z(d_pz);  
  Enthalpy t1(t_x, g_x, u_x, time, NULL, &f_x);
  Enthalpy t2(t_y, g_y, u_y, time, NULL, &f_y);
  Enthalpy t3(t_z, g_z, u_z, time, NULL, &f_z);

  /*-----------------------+
  |  discretize and solve  |
  +-----------------------*/
  t1.discretize();    t1.direct();
  t2.discretize();    t2.direct();
  t3.discretize();    t3.direct();

  boil::plot->plot(t_x,"t1",0);
  boil::plot->plot(t_y,"t2",0);
  boil::plot->plot(t_z,"t3",0);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(t_z, "test", 0);
}	
