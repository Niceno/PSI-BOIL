#include "Include/psi-boil.h"
#include <vector>
#include "update_step.cpp"
#define PHASECHANGE

/* parameters */
const int NX = 64;
const int NY = NX;
const int NZ = 4*NX;
const int NINJ = 12;

const real LX = 0.02;
const real LY = LX/real(NX)*real(NY);
const real LZ = LX/real(NX)*real(NZ);

/* constants */
const real gravity =9.8;
const real Tbulk = 20.0;
const real Tsat  = 100.0;

/****************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  if(argc==1){
    boil::oout<<"One command line argument is required!"<<"\n";
    boil::oout<<"./Boil wmin (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*--------+
  |  grids  |
  +--------*/
  Grid1D gx(Range<real>( -LX/2.0, LX/2.0 ), NX, Periodic::no());
  Grid1D gy(Range<real>( -LY/2.0, LY/2.0 ), NY, Periodic::no());
  Grid1D gz(Range<real>( 0.0,  LZ ), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  int iref = d.local_i(NX/2);
  int jref = d.local_j(NY/2);
  int kref = d.local_k(1);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar press(d), p  (d), f  (d); // pressure
  Scalar c  (d), g  (d), kappa(d), cold(d); // concentration
  Scalar tpr(d), q  (d);           // temperature
  Scalar step(d), sflag(d);        // step
  Scalar mdot(d);                  // phase-change rate
  Scalar mu_t(d);                  // eddy viscosity
  Scalar rgnid(d);  		   //to track regions

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd(Range<int>(NX/2+1-NINJ/2,NX/2+NINJ/2)
                        , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                        , Dir::kmin()
                        , BndType::inlet(),0.0,0.0,1.0 ) );
    uvw.bc(m).add( BndCnd(Range<int>(1,NX)
                        , Range<int>(1,NY/2+1-NINJ/2-1)
                        , Dir::kmin()
                        , BndType::wall() ) );
    uvw.bc(m).add( BndCnd(Range<int>(1,NX/2+1-NINJ/2-1)
                        , Range<int>(NY/2-NINJ/2,NY/2+NINJ/2)
                        , Dir::kmin()
                        , BndType::wall() ) );
    uvw.bc(m).add( BndCnd(Range<int>(NX/2+NINJ/2+1,NX)
                        , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                        , Dir::kmin()
                        , BndType::wall() ) ); 
    uvw.bc(m).add( BndCnd(Range<int>(1,NX)
                        , Range<int>(NY/2+NINJ/2+1,NY)
                        , Dir::kmin()
                        , BndType::wall() ) ); 

    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
  }

  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press = p.shape();
  mdot = p.shape();
  mu_t = p.shape();
  f = p.shape();
  g = p.shape();

  c.bc().add( BndCnd(Range<int>(NX/2+1-NINJ/2,NX/2+NINJ/2)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::dirichlet(), 0.0 ) );
  c.bc().add( BndCnd(Range<int>(1,NX)
                   , Range<int>(1,NY/2+1-NINJ/2-1)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd(Range<int>(1,NX/2+1-NINJ/2-1)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd(Range<int>(NX/2+NINJ/2+1,NX)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd(Range<int>(1,NX)
                   , Range<int>(NY/2+NINJ/2+1,NY)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  step = c.shape();
  sflag = c.shape();
  rgnid = c.shape();

  tpr.bc().add( BndCnd(Range<int>(NX/2+1-NINJ/2,NX/2+NINJ/2)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::dirichlet(), 151.8 ) );
  tpr.bc().add( BndCnd(Range<int>(1,NX)
                   , Range<int>(1,NY/2+1-NINJ/2-1)
                   , Dir::kmin()
                   , BndType::neumann() ) );
  tpr.bc().add( BndCnd(Range<int>(1,NX/2+1-NINJ/2-1)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::neumann() ) );
  tpr.bc().add( BndCnd(Range<int>(NX/2+NINJ/2+1,NX)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::neumann() ) );
  tpr.bc().add( BndCnd(Range<int>(1,NX)
                   , Range<int>(NY/2+NINJ/2+1,NY)
                   , Dir::kmin()
                   , BndType::neumann() ) );

  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), Tbulk ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), Tbulk ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), Tbulk ) );
  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), Tbulk ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), Tbulk ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d);
  vapor  .mu    (1.396e-5);
  vapor  .rho   (2.668);
  vapor  .cp    (2410*2.668);
  vapor  .lambda(0.03025);
  liquid.mu    (0.001);
  liquid.rho   (1000.0);
  liquid.cp    (4184.0*1000.0);
  liquid.lambda(0.600);

  Matter mixed(liquid, vapor, & step);
  mixed.sigma(5.9e-2);
  const real latent=2258.0*1e3;
  const real liquid_drhodt=-0.7;  //[kg/m3K]
  const real vapor_drhodt=-0.0017; //[kg/m3K]

  /*-------+
  |  Time  |
  +-------*/
  const int ndt = 50000000;  // total time step to be computed
  const real tint = 1.0e-3;  // time interval of tecplot file
  const int nint  = 5000;    // time interval of *.bck file
  const real dxmin = d.dxyz_min();
  const real dt  = 1.0*pow(vapor.rho()->value()*pow(dxmin,3.0)
                   / (2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const real cfl_limit = 0.2;
  Times time(ndt, dt);
  time.print_time(false);
  time.set_coef_dec(0.2);

 /*----------------------------+
  |  Tracking algorithm object |
  +---------------------------*/
  Floodfill idrgns(c, rgnid, &uvw, time);
  idrgns.set_out_freq(10);
  boil::oout<<"main:Floodfill:output frequency= "<<idrgns.get_out_freq()<<"\n";

  /*---------+
  |  solver  |
  +---------*/
  /* linear solver */
  Krylov * solver = new CG(d, Prec::ic2());
  /* pressure */
  Pressure pr( p,   f,   uvw, time, solver, &mixed );
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);
  /* momentum */
  Momentum ns( uvw, xyz,      time, solver, &mixed );
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::upwind()); 
  /* color function */
  CIPCSL2 conc(c, g, kappa, uvw, time, solver);
  conc.set_itsharpen(10);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_cangle(90.0);  // hydrophilic
  /* enthalpy */
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, & mixed ,Tsat);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  /* phase change */
  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw, time, &mixed, latent, Tsat);

  /*eddy viscosity*/
  Model mod;

  /*---------------------------------------------------------------------------+
  ! start computation                                                          |
  +---------------------------------------------------------------------------*/
  /*-------------------+
  |  check if restart  |
  +-------------------*/
  std::fstream input;
  int irun = 0;
  if(boil::cart.iam()==0){
    input.open("run.txt", std::ios::in);
    if( !input.fail() ) {
      input >> irun;
      std::cout<<"read irun.  irun= "<<irun<<"\n";
    }
    input.close();
  }
  boil::cart.sum_int(&irun);
  if (irun==1){
    boil::oout<<"exit job due to irun=1"<<"\n";
    exit(0);
  }

  if(boil::cart.iam()==0){
    std::fstream output;
    output.open("run.txt", std::ios::out);
    output << 1 << boil::endl;
    output.close();
  }

  int ts=0;
  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    uvw.  load("uvw",ts);
    press.load("press",ts);
    conc. load("conc",ts);
    tpr.  load("tpr",ts);
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
  } else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    c = 1.0;
    c.bnd_update();
    c.exchange_all();
    conc.init();

    tpr = Tbulk;
    tpr.bnd_update();
    tpr.exchange_all();

    boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,rgnid
                   ,"uvw-c-tpr-press-mdot-mu_t-rgnid",0);
  }
  input.close();

  update_step(c, step, sflag);
  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  Time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "##################" << boil::endl;

    /*------------------+
    |  reset body force |
    +------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /*-----------------+
    |  eddy viscosity  |
    +-----------------*/
    mod.smagorinsky( &ns, &mu_t, 0.173 );

    /*---------------+
    |  phase change  |
    +---------------*/
#ifdef PHASECHANGE
    pc.update(&mu_t);
    ns.vol_phase_change(&f);
#endif

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* gravity force */
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k){
      real phil=step[i][j][k];
      real phiv=1.0-phil;
      real deltmp=tpr[i][j][k]-Tsat;
      real rhomix = (liquid.rho()->value() + liquid_drhodt*deltmp)*phil
                  + (vapor.rho()->value()  + vapor_drhodt*deltmp)*phiv;
      if(d.ibody().on(m,i,j,k))
        xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
    }
    xyz.exchange();

    /* surface tension */
    conc.tension(&xyz, mixed, step);

    /* eddy viscosity near outlet to avoid divergence */
    const real z0 = 0.9*LZ;
    const real z1 = 0.95*LZ;
    for_avijk(c,i,j,k){
      real ztmp = fabs(c.zc(k));
      if(ztmp>z0){
        real coef = std::min((ztmp-z0)/(z1-z0),1.0);
        mu_t[i][j][k] += coef * liquid.mu()->value() * 100;
      }
    }

    ns.discretize( &mu_t );
    pr.discretize();
    pr.coarsen();

    /* intermediate velocity */
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-14));

    /* pressure */
    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-6))) OMS(converged);
    p.exchange();

    /* update velocity and pressure */
    ns.project(p);
    press += p;

    /* shift pressure */
    real pref = 0.0;

    if(kref>0){
      if(jref>0){
        if(iref>0){
          pref = press[iref][jref][kref];
        }
       }
     }

    boil::cart.sum_real(&pref);

    for_vijk(press,i,j,k){
        press[i][j][k] -= pref;
    }
    press.bnd_update();
    press.exchange_all();

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    cold = c;
    conc.advance();
    pc.modify_vel(uvw,c,cold);

    /* outlet region: delete bubbles artificially */
    const real z0b = 0.90*LZ;
    const real z1b = 0.95*LZ;
    for_avk(c,k){
      if(c.zc(k)>z0b){
        real coef=std::min((c.zc(k)-z0b)/(z1b-z0b),1.0);
        for_avij(c,i,j){
          c[i][j][k]= (1.0-coef)*c[i][j][k] + coef*1.0;
          tpr[i][j][k]= Tbulk;
        }
      }
    }
    c.bnd_update();
    c.exchange_all();
    conc.update_node(c);
    update_step(c, step, sflag);
    tpr.bnd_update();
    tpr.exchange_all();


    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
#ifdef PHASECHANGE
    enthFD.discretize( &mu_t );
    enthFD.new_time_step( &mu_t );
    enthFD.solve(ResRat(1e-16),"enthFD");
#endif


    /*------------------------------------------+
    |  Identify region with tracking algorithm  |
    +------------------------------------------*/
    idrgns.identify_regions();

    /*-------------+
    |  dt control  |
    +-------------*/
    real cflmax = ns.cfl_max();
    time.control_dt(cflmax, cfl_limit, dt);

    if (time.dt()<1.0e-9) {
      boil::oout<<"Too small dt: "<<dt<<"\n";
      exit(0);
    }

    /*--------------+
    |  output data  |
    +--------------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      tpr.exchange_all();
      boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,rgnid
                     ,"uvw-c-tpr-press-mdot-mu_t-rgnid",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if((time.current_step()) % (nint)==0 ) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
      if( boil::cart.iam()==0) {
        std::fstream output;
        std::stringstream ss;
        ss <<"time-"<<time.current_step()<<".txt";
        std::string fname = ss.str();
        int len = fname.length();
        char * cfname = new char[len+1];
        memcpy(cfname, fname.c_str(), len+1);
        output << std::setprecision(16);
        output.open(cfname, std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output.close();
      }
    }

    if( boil::timer.current_min() > (wmin-30.0)
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
      std::fstream output;
      output << std::setprecision(16);
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      output.open("run.txt", std::ios::out);
      output << 0 << boil::endl;
      output.close();
      boil::timer.stop();
      boil::timer.report();
      uvw  .rm("uvw", ts);
      press.rm("press", ts);
      conc .rm("conc", ts);
      tpr  .rm("tpr", ts);
      exit(0);
    }
  }

  boil::timer.stop();
  boil::timer.report();
}
