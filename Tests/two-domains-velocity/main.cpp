#include "Include/psi-boil.h"
#include <iomanip>
#include "update_step.cpp"

/* boundary conditions */
const real LX =   0.1;
const real LY =   0.02;
const int  gLevel = 2;
const int  NX =  32*gLevel;
const int  NY =  16*gLevel;
const int  NZ =  32*gLevel;

const real tsat = 100.0;
const real gravity = 9.8;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  if(argc==1){
    boil::oout<<"One command line argument is required!"<<"\n";
    boil::oout<<"./Boil wmin (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  real       dpdx  =-1.0e+2;

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  //boil::plot = new PlotTEC(AsNodes::no(),Buffers::yes());
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D g_x1( Range<real>(0, LX), 
               NX, Periodic::yes());
  Grid1D g_x2( Range<real>(LX, 2.0*LX), 
               NX, Periodic::no());

  Grid1D g_y( Range<real>(-0.5*LY, 0.5*LY),  
              NY, Periodic::yes());

  Grid1D g_z( Range<real>(-0.005,LY-0.005),  
              Range<real>( LY/real(NZ), LY/real(NZ) ),
              NZ, Periodic::no());

  /*----------+
  |  domains  |
  +----------*/
  Body floor_1("floor.stl");
  Body floor_2("floor.stl");
  Domain dom_1(g_x1, g_y, g_z, & floor_1);
  Domain dom_2(g_x2, g_y, g_z, & floor_2);

  /*-----------------------------------+
  |  define unknowns and default b.c.  |
  -----------------------------------*/
  Vector uvw_1(dom_1), xyz_1(dom_1); // vel 1
  Vector uvw_2(dom_2), xyz_2(dom_2); // vel 2
  Scalar press_1(dom_1), p_1  (dom_1), f_1  (dom_1); // p 1
  Scalar press_2(dom_2), p_2  (dom_2), f_2  (dom_2); // p 2
  Scalar t_1  (dom_1), q_1  (dom_1); // t 1
  Scalar t_2  (dom_2), q_2  (dom_2); // t 2
  Scalar c_1  (dom_1), g_1  (dom_1), step_1(dom_1), sflag(dom_1); // c 1
  Scalar c_2  (dom_2), g_2  (dom_2), step_2(dom_2);               // c 2
  Scalar wd_1 (dom_1), ws_1 (dom_1), mu_t_1(dom_1); // wall distance 1
  Scalar wd_2 (dom_2), ws_2 (dom_2), mu_t_2(dom_2); // wall distance 2
  Scalar mdot_2  (dom_2); // phase change mdot 2

  real ** copy_plane;
  real *** copy_planeVec;

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  #include "boundary_1.cpp"
  #include "boundary_2.cpp"

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  #include "property_1.cpp"
  #include "property_2.cpp"

  const real latent=2258.0*1e3;
  const real liquid_drhodt=-0.7;   //[kg/m3K]
  const real vapor_drhodt=-0.0017; //[kg/m3K]

  /*-------------------+
  |  time integration  |
  +-------------------*/
  const real dxmin = std::max(dom_1.dxyz_min(), dom_2.dxyz_min());
  const real dt  =10.0*pow(vapor_1.rho()->value()*pow(dxmin,3.0)
                 / (2.0*3.1415*mixed_1.sigma()->value()),0.5);
  const real cfl_limit=0.25;
  Times time(100000, dt); /* ndt, dt */
  const real tint = 1.0e-2;
  int nint = 5000;
  time.print_time(false);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  #include "solver_1.cpp"
  #include "solver_2.cpp"

  /*----------+
  |  restart  |
  +----------*/
  bool restart=false;
  int ts=0;
  #include "restart.cpp"
  #include "restart_1.cpp"
  #include "restart_2.cpp"

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  #include "initial_1.cpp"
  #include "initial_2.cpp"

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


    /*------------------------------------------------------------------------+
    |  domain 1                                                               |
    +------------------------------------------------------------------------*/
    /*---------------------+
    |  reset source terms  |
    +---------------------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz_1,m,i,j,k)
        xyz_1[m][i][j][k] = 0.0;
    /* heat source */
    q_1=0.0;

    /*-----------+
    |  momentum  |
    +-----------*/
    /* gravity force */
    Comp m = Comp::w();
    for_vmijk(xyz_1,m,i,j,k){
      real phil=step_1[i][j][k];
      real phiv=1.0-phil;
      real deltmp=t_1[i][j][k]-tsat;
      real rhomix = (liquid_1.rho()->value() + liquid_drhodt*deltmp)*phil
                  + (vapor_1.rho()->value()  + vapor_drhodt*deltmp)*phiv;
      if(dom_1.ibody().on(m,i,j,k))
        xyz_1[m][i][j][k] += -gravity * xyz_1.dV(m,i,j,k) * rhomix;
    }
    /* pressure gradient */
    m = Comp::u();
    for_vmijk(xyz_1,m,i,j,k){
      if(dom_1.ibody().on(m,i,j,k))
        xyz_1[m][i][j][k] += -dpdx * xyz_1.dV(m,i,j,k);
    }

    /*surface tension */
    conc_1.tension(&xyz_1, mixed_1, step_1);

    /* turbulence model */
    tm_1.wale( & ns_1, & mu_t_1, 0.325 );
    tm_1.tau_wall( & ns_1, wd_1, & xyz_1 );

    /* intermediate velocity */
    ns_1.discretize( &mu_t_1 );
    pr_1.discretize();
    pr_1.coarsen();
    ns_1.new_time_step();
    ns_1.grad(press_1);
    ns_1.solve(ResRat(1.0e-4));

    /* pressure poisson equation */
    p_1 = 0.0;
    sol_1.vcycle(ResRat(1.0e-4));

    /* update velocity and pressure */
    ns_1.project(p_1);
    press_1 += p_1;
    press_1.bnd_update();
    press_1.exchange_all();

    /* shift pressure */
    real pmin_1=1.0e+300;
    for_vijk(press_1,i,j,k){
      if(dom_1.ibody().on(i,j,k)){
        if(pmin_1>press_1[i][j][k]) pmin_1=press_1[i][j][k];
      }
    }
    boil::cart.min_real(&pmin_1);

    for_vijk(press_1,i,j,k){
      if(dom_1.ibody().on(i,j,k)){
        press_1[i][j][k] -= pmin_1;
      } else {
        press_1[i][j][k] = 0.0;
      }
    }
    press_1.bnd_update();
    press_1.exchange_all();

    /*-----------------+
    |  color function  |
    +-----------------*/
    conc_1.advance();
    update_step(c_1, step_1, sflag);

    /*--------------+
    |  temperature  |
    +--------------*/
    en_1.discretize( &mu_t_1 );
    en_1.new_time_step( &mu_t_1 );
    en_1.solve(ResRat(1.0e-12));

    /*------------------------------------------------------------------------+
    |  domain 2                                                               |
    +------------------------------------------------------------------------*/
    /*---------------------+
    |  copy domain 1 to 2  |
    +-------------------- */
    uvw_1.bnd_extract( Dir::imax(), &copy_planeVec );
    uvw_2.bnd_insert ( Dir::imin(),  copy_planeVec );
    t_1.bnd_extract( Dir::imax(), &copy_plane );
    t_2.bnd_insert ( Dir::imin(),  copy_plane );
    conc_1.bnd_extract( Dir::imax(), &copy_plane );
    conc_2.bnd_insert ( Dir::imin(),  copy_plane );
    update_step(c_2, step_2, sflag);

    /*---------------------+
    |  reset source terms  |
    +---------------------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz_2,m,i,j,k)
        xyz_2[m][i][j][k] = 0.0;
    /* heat source */
    q_2=0.0;

    /* turbulence model */
    tm_2.wale( & ns_2, & mu_t_2, 0.325 );
    tm_2.tau_wall( & ns_2, wd_2, & xyz_2 );

    /*---------------+
    |  phase change  |
    +---------------*/
    pc_2.update(&mu_t_2);
    //pc.micro(&xyz_2);
    //update_step(c, step, sflag);  // 0119 need to update step for with IB
    ns_2.vol_phase_change(&f_2);

    /*-----------+
    |  momentum  |
    +-----------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz_2,m,i,j,k)
        xyz_2[m][i][j][k] = 0.0;

    /* gravity force */
    m = Comp::w();
    for_vmijk(xyz_2,m,i,j,k){
      real phil=step_2[i][j][k];
      real phiv=1.0-phil;
      real deltmp=t_2[i][j][k]-tsat;
      real rhomix = (liquid_2.rho()->value() + liquid_drhodt*deltmp)*phil
                  + (vapor_2.rho()->value()  + vapor_drhodt*deltmp)*phiv;
      if(dom_2.ibody().on(m,i,j,k))
        xyz_2[m][i][j][k] += -gravity * xyz_2.dV(m,i,j,k) * rhomix;
    }

    /* pressure gradient */
    m = Comp::u();
    for_vmijk(xyz_2,m,i,j,k){
      if(dom_2.ibody().on(m,i,j,k))
        //xyz_2[m][i][j][k] += -dpdx * xyz_2.dV(m,i,j,k);
        xyz_2[m][i][j][k] += 0.0;
    }

    /* surface tension */
    conc_2.tension(&xyz_2, mixed_2, step_2);

    /* intermediate velocity */
    ns_2.discretize( &mu_t_2 );
    pr_2.discretize();
    pr_2.coarsen();
    ns_2.new_time_step();
    ns_2.grad(press_2);
    ns_2.solve(ResRat(1.0e-4));

    /* pressure poisson equation */
    p_2 = 0.0;
    sol_2.vcycle(ResRat(1.0e-4));

    /* update velocity and pressure */
    ns_2.project(p_2);
    press_2 += p_2;
    press_2.bnd_update();
    press_2.exchange_all();

    /* shift pressure */
    real pmin_2=1.0e+300;
    for_vijk(press_2,i,j,k){
      if(dom_2.ibody().on(i,j,k)){
        if(pmin_2>press_2[i][j][k]) pmin_2=press_2[i][j][k];
      }
    }
    boil::cart.min_real(&pmin_2);

    for_vijk(press_2,i,j,k){
      if(dom_2.ibody().on(i,j,k)){
        press_2[i][j][k] -= pmin_2;
      } else {
        press_2[i][j][k] = 0.0;
      }
    }
    press_2.bnd_update();
    press_2.exchange_all();

    /*-----------------+
    |  color function  |
    +-----------------*/
    conc_2.advance();
    update_step(c_2, step_2, sflag);

    /*--------------+
    |  temperature  |
    +--------------*/
    en_2.discretize( &mu_t_2 );
    en_2.new_time_step( &mu_t_2 );
    en_2.solve(ResRat(1.0e-12));

    /*-------------------------------------------------------------------------+
    |  dt control                                                              |
    +-------------------------------------------------------------------------*/
    real cflmax_1 = ns_1.cfl_max();
    real cflmax_2 = ns_2.cfl_max();
    real cflmax = std::max(cflmax_1, cflmax_2);
    time.control_dt(cflmax, cfl_limit, dt);
    boil::oout<<"main:cflmax= "<<cflmax_1<<" "<<cflmax_2<<"\n";

    /*---------+
    |  output  |
    +---------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw_1,c_1,t_1,press_1,mu_t_1,
                     "d1-uvw-c-t-press-mu_t", iint);
      boil::plot->plot(uvw_2,c_2,t_2,press_2,mu_t_2,mdot_2,
                     "d2-uvw-c-t-press-mu_t-mdot", iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /*---------+
    |  backup  |
    +---------*/
    #include "backup_1.cpp"
    #include "backup_2.cpp"
    #include "backup.cpp"

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}
