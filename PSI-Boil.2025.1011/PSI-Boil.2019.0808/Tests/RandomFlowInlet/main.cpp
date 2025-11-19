#include "Include/psi-boil.h"
#include <iomanip>

/* average velocity */
const real vel_ave = 20.0;   // (m/s)

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*---------------------+
  |  output file format  |
  +---------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  const real Z0 = -0.56e-3;
  const real Z1 =  6.475e-3;    //       Z0 <= z <= Z1
  const real LX =  Z1;          //        0 <= x <= Z1
  const real LY =  Z1*0.5;      // -0.25*Z1 <= y <= 0.25*Z1
  const int gLevel = 1;
  const int NX  = 64 * gLevel;
  const int NY  = NX * LY / LX;
  const int NZ0 =  4 * gLevel;
  const int NZ1 = 60 * gLevel;
  const real dx = LX/real(NX);

  Grid1D gx( Range<real>(0, LX), NX, Periodic::no());
  Grid1D gx_in(Range<real>(0, dx*2.0 ), 2, Periodic::no());

  Grid1D gy( Range<real>(-0.5*LY, 0.5*LY), NY, Periodic::yes());

  Grid1D gz0(Range<real>( Z0, 0 ), NZ0, Periodic::no());
  Grid1D gz1(Range<real>( 0, Z1 )
           , Range<real>( 0.5*dx, 2.0*dx)
           , NZ1, Periodic::no());
  Grid1D gz (gz0, gz1, Periodic::no());

  /*----------+
  |  domains  |
  +----------*/
  Body floor("floor.stl");
  Domain dom(gx, gy, gz, & floor);
  Domain d_in(gx_in, gy, gz, "inlet_dom", Decompose::no());

  /*-----------------------------------+
  |  define unknowns and default b.c.  |
  -----------------------------------*/
  Vector uvw(dom), xyz(dom); // vel 1
  Scalar press(dom), p(dom), f(dom); // p 1
  Scalar wd(dom), ws(dom), mu_t(dom); // wall distance 1

  Vector uvw_in(d_in);   // inlet velocity
  real uvw_bulk[3] = {vel_ave, 0.0, 0.0};
  real ** copy_plane;
  real *** copy_planeVec;

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::insert() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  press.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 0.0 ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  p.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 0.0 ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  f    = p.shape();
  mu_t = p.shape();

  mu_t.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  mu_t.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  mu_t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  mu_t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  mu_t.bc().add( BndCnd( Dir::kmin(), BndType::neumann()  ) );
  mu_t.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );
  wd = mu_t.shape();

  for_m(m) {
    uvw_in.bc(m).add( BndCnd( Dir::imin(), BndType::insert() ) );
    uvw_in.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw_in.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw_in.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw_in.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw_in.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(dom);
  vapor  .mu    (0.129e-04);
  vapor  .rho   (1.076e+00);
  vapor  .cp    (2112.3*1.076e+00);
  vapor  .lambda(0.268e-01);

  /*-------------------+
  |  time integration  |
  +-------------------*/
  const real dxmin = dom.dxyz_min();
  const real dt  = 1.0e-6;
  const int ndt = 100;
  Times time(ndt, dt); /* ndt, dt */
  const real tint = 1.0e-5;
  const real cfl_limit=0.25;
  time.set_coef_dec(0.25);
  time.print_time(false);

  /*-----------------+
  |  define solvers  |
  +-----------------*/

  Krylov * solver = new CG(dom, Prec::di());

  /* momentum */
  Momentum ns( uvw, xyz, time, solver, &vapor );
  ns.convection_set(TimeScheme::forward_euler());

  /* pressure poisson equation */
  Pressure pr( p, f, uvw, time, solver, &vapor );
  AC sol( &pr );
  sol.stop_if_diverging(true);
  sol.min_cycles(3);
  sol.max_cycles(8);

  /* wall distance */
  Distance di(wd,  ws,  uvw, time, solver);
  di.compute();
  boil::plot->plot(wd, "wd");

  /* eddy viscosity */
  Model tm;

  /* random flow generator */
  int ne = 100;
  real tl = 1e-4; // turbulent length scale :original 0.1
  real tt = 1e-4; // turbulent time scale   :original 0.1
  RandomFlow rf( ne, tt, tl );

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  Comp m = Comp::u();
  for_avmijk(uvw,m,i,j,k){
    if (uvw.zc(m,k)<0.0){
      uvw[m][i][j][k]=0.0;
    } else {
      uvw[m][i][j][k]=vel_ave;
    }
  }

   boil::plot->plot(uvw,press,mu_t, "uvw-press-mu_t", 0);

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# DT:        " << time.dt() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "########################" << boil::endl;

    /*---------------------+
    |  reset source terms  |
    +---------------------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /*------------------+
    |  turbulence model |
    +------------------*/
    tm.smagorinsky( & ns, & mu_t, 0.173);

    /*-----------+
    |  momentum  |
    +-----------*/
    /* inlet plane */
    uvw_in.randomize(rf, tl/tt, time.current_time());
    uvw_in += uvw_bulk;
    for_m(m){
      for_vmijk(uvw_in,m,i,j,k){
        real zz = uvw_in.zc(m,k);
        if(zz<0.0) {
          uvw_in[m][i][j][k] = 0.0;
        }
      }
    }
    uvw_in.bnd_extract_global( Dir::imin(), &copy_planeVec, &uvw );
    uvw.bnd_insert ( Dir::imin(),  copy_planeVec );

    ns.discretize( &mu_t );
    /* wall function */
    tm.tau_wall( & ns, wd, & xyz );
    pr.discretize();
    pr.coarsen();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1.0e-8));

    p = 0.0;
    sol.vcycle(ResRat(1.0e-6));
    ns.project(p);
    press += p;
    press.bnd_update();
    press.exchange_all();

    /*-------------------------------------------------------------------------+
    |  dt control                                                              |
    +-------------------------------------------------------------------------*/
    real cflmax = ns.cfl_max();
    boil::oout<<"main:cflmax= "<<cflmax<<"\n";
    time.control_dt(cflmax, cfl_limit, dt);

    /*---------+
    |  output  |
    +---------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,press,mu_t, "uvw-press-mu_t", iint);
      iint = int(time.current_time()/tint) + 1;
    }

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}
