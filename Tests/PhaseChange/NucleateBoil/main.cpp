#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#include "update_step.cpp"
//#define USE_VOF   // If USE_VOF is defined, then VOF is used.
                    // Otherwise, CIPCSL2 is used for interface tracking method.

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Simplified version of
* Sato, Niceno, JCP 300 (2015) 20-52
* Duan's case: Case 1 in Table 3
* Grid: Grid e (gLevel=3) in Table 4
* The result may differ from the paper due to the simplification
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* domain dimensions */
const real LX1 = 0.003;
const real LX2 = 0.0075;
const real LZ0 = -(250.0)*1.0e-6;
const real LZ1  = 0.008;
const real LZ2  = 0.020;
const int gLevel = 3;  //grid level=2,3,4
const int NX1  = 12*gLevel;
const int NX2  = 10*gLevel;
const int NZ0  =  4*gLevel;  // solid
const int NZ1  = 32*gLevel;
const int NZ2  = 20*gLevel;

/* parameter for boundary or initial condition */
const real tsat   = 1.0E+02;   // saturation temperature [deg.]
const real tout  = tsat - 0.5; // temperature at outlet and side boundaries [deg.]
const real tact  = tsat + 9.0; // nucleation activation temperature [deg.]

/* constants */
const real gravity = 9.8;

/* heater */
const real h_width = 0.010;       // heater width
const real qflux   = 28.7*1000.0; // heater power[W/m2]

/******************************************************************************/
int main(int argc, char ** argv) {

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

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx0( Range<real>(-LX2,-LX1),
              Range<real>(2.0*LX1/real(NX1),LX1/real(NX1)),
              NX2, Periodic::no() );
  Grid1D gx1( Range<real>(-LX1,LX1), 2*NX1, Periodic::no() );
  Grid1D gx2( Range<real>( LX1, LX2),
              Range<real>(LX1/real(NX1),2.0*LX1/real(NX1)),
              NX2, Periodic::no() );
  Grid1D gt1 (gx0 , gx1, Periodic::no());
  Grid1D gx  (gt1 , gx2, Periodic::no());
  Grid1D gz0( Range<real>(LZ0, 0.0), NZ0, Periodic::no() );
  Grid1D gz1( Range<real>(0.0, LZ1), NZ1, Periodic::no() );
  Grid1D gz2( Range<real>(LZ1, LZ2),
              Range<real>(LZ1/real(NZ1),4.0*LZ1/real(NZ1)),
              NZ2, Periodic::no() );
  Grid1D gztmp ( gz0, gz1, Periodic::no());
  Grid1D gz    ( gztmp, gz2, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Domain d(gx, gx, gz, & floor);  // immerse boundary method

  /*------------------+
  |  define unknowns  |
  +------------------*/
  // Note: character in "" will appear in TEC file, e.g. "vfl"
  Vector uvw(d), xyz(d);               // uvw: velocity [m/s], xyz: source term for NS (body force) [N]
  Scalar p(d), f(d), press(d,"press"); // p: pressure difference in time step [Pa]
				       // f: source term of pressure Poisson equation [m3/s2]
  Scalar c(d,"vfl"), g(d), kappa(d,"kappa"); // c: color function (= volume fraction of liquid) [-]
					     // g: source term of the equation for color [1/s]
					     // kappa: curvature [1/m]
  Scalar tpr(d,"temperature"), q  (d); // tpr: temperature [deg.], q: heat source [W]
  Scalar mdot(d,"mdot");               // mdot: phase-change rate [kg/m3s]
  Scalar dmicro(d,"dmicro");           // dmicro: micro-layer film thickness [m]
  Scalar mu_t(d,"mut");                // eddy viscosity [Pa.s]
#ifndef USE_VOF
  Scalar step(d), sflag(d);            // step: step function for color function [-]
#endif

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /* copy b.c. from p */
  press = p.shape();
  f     = p.shape();
  mdot  = p.shape();
  q     = p.shape();
  mu_t  = p.shape();
  kappa = p.shape();
  g     = p.shape();
  dmicro = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
#ifndef USE_VOF
  step = c.shape();
  sflag = c.shape();
#endif

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), tout ) );

  /* heater power */
  real qsrc;
  for_vk(tpr,k){
    if(approx(tpr.zn(k+1),0.0)){
      qsrc=qflux/tpr.dzc(k);  // [W/m3]
    }
  }
  boil::oout<<"#qsrc= "<<qsrc<<"\n";

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d), sapphire(d);
  vapor  .mu    (1.255e-5);  // dynamic viscosity [Pa.s]
  vapor  .rho   (0.597);     // density [kg/m3]
  vapor  .cp    (2030*0.597);// heat capacity [J/m3K]
  vapor  .lambda(0.025);     // thermal conductivity [W/mK]
  liquid.mu    (0.28e-3);
  liquid.rho   (958.4);
  liquid.cp    (4215.9*958.4);
  liquid.lambda(0.679);
  sapphire.rho    (3980.0);
  sapphire.cp     (750.0*3980.0);
  sapphire.lambda (35.0);
  const real latent=2258.0*1e3; // latent heat [J/kg]
#ifdef USE_VOF
  Matter mixed(liquid, vapor, & c);  // property averaged with c
#else
  Matter mixed(liquid, vapor, & step);  // property averaged with step
#endif
  mixed.sigma(5.9e-2);  // surface tension coefficient [N/m]
  const real liquid_drhodt=-0.7;   // temperature derivation of liquid density [kg/m3K]
  const real vapor_drhodt=-0.0017; // temperature derivation of vapor density [kg/m3K]

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 500000;   // total number of time step to be computed
  const real tint = 5.0e-4;  // time interval for TEC file [s]
  const int  nint = 1000*gLevel;  // time step interval for backup file
  const real dxmin = d.dxyz_min();// minimum cell size in fluid domain
  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<"\n";
  const real dt  =10.0*pow(0.5*pow(dxmin,3.0)  // 10 is acceleration
                 / (2.0*3.1415*mixed.sigma()->value()),0.5);
  const real cfl_with = 0.05;  // max CFL number when liquid-vapor interface exists
  const real cfl_wo   = 0.2;   // max CFL number when no liquid-vapor interface exists
  Times time(ndt, 0.002);      // 0.002 is initial time increment
  time.print_time(false);
  time.set_coef_dec(0.75);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::ic2());
  Krylov * solver2 = new CG(d, Prec::di());

  /* color function */
#ifndef USE_VOF
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_itsharpen(10);
#else
  VOF conc  (c,  g, kappa, uvw, time, solver);
#endif
  conc.set_cangle(0.0);  // static contact angle

  /* momentum equation */
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.convection_set(TimeScheme::forward_euler());

  /* pressure solver */
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);  // stop iteration when residual increases
  multigrid.min_cycles(3);            // min no. iteration for multigrid
  multigrid.max_cycles(10);           // max no. iteration for multigrid

  /* enthalpy equation */
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver2, & mixed ,tsat, & sapphire);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  /* nucleation site */
  real rseed = 1.5 * dxmin; // radius of seed bubble [m]
  real zplant=-1.0; // when bottom of bubble reaches zplant, next seed is set
  const real dmicro_min=1.0e-8;  // min microlayer thickness to avoid divided-by-zero
  Nucleation nucl( &c, &tpr, &q, &time, dmicro, &mixed, rseed,
                   dmicro_min, latent, conc.get_cangle(), tsat, &sapphire);
  nucl.set_slope(1.0*4.46e-3);  // C_slope, initial microlayer thickness = C_slope * r
  nucl.set_seed_period(1e-5);   // bubble is enforced to be planted during this period
				// after activation of the site [s]

  /* phase change */
#ifdef USE_VOF
  PhaseChange pc( mdot, tpr, q, c, g, f, c, uvw, time, &mixed , latent,
                  tsat, &sapphire, &nucl);
#else
  PhaseChange pc( mdot, tpr, q, c, g, f, step, uvw, time, &mixed , latent,
                  tsat, &sapphire, &nucl);
#endif

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
    uvw.  load("uvw",   ts);
    press.load("press",   ts);
#ifndef USE_VOF
    conc. load("conc", ts);
#else
    c.load("c", ts);
#endif
    tpr.  load("tpr", ts);
    nucl. load("nucl", ts);
    dmicro.load("dmicro", ts);
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
    /* set allow_replant */
    if( ns.cfl_max()<=cfl_with && time.dt()<=dt) {
      boil::oout<<"Restart: Allow replant.\n";
      for(int nsd=0; nsd<nucl.size(); nsd++){
        nucl.sites[nsd].set_allow_replant(true);
      }
    } else {
      boil::oout<<"Restart: Deny replant for this step.\n";
    }
  } else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    /*--------------------+
    |  initial condition  |
    +--------------------*/
    // make domain full of liquid
    for_vijk(c,i,j,k) 
      c[i][j][k] = 1.0;

    c.bnd_update();   // after Scalar or Vector variable is enforced to be modified in main.cpp,
    c.exchange_all(); // bnd_update() and exchnage_all()/exchange() must be called.

    // initialization of VOF/CIPCSL
    conc.init();

    /* add nucleation site */
    real zsite=rseed*cos(45.0/180.0*boil::pi);
    nucl.add(Site( 0.000,  0.000, zsite, tact, zplant)); // center of spherical seed bubble is
							 // (0.0, 0.0, zsite)
							 // if zsite = 0, hemi-sphere.

#if 1
    tpr = tout;
#else
    // temperature with thermal boundary layer
    const real twall = tsat + 8.95;// used for initial temperature field [deg.]
    const real ztconst = 1.0*1.0e-3;  // thermal boundary layer thickness
    for_vijk(c,i,j,k) {
      if(tpr.zc(k)<0.0) {
        tpr[i][j][k] = twall;
      } else if (tpr.zc(k)<ztconst) {
        tpr[i][j][k] = twall + (tout-twall)/ztconst * tpr.zc(k);
      } else {
        tpr[i][j][k] = tout;
      }
    }
#endif
    tpr.bnd_update();
    tpr.exchange_all();

    // output of initial condition
    boil::plot->plot(uvw,c,tpr,press,mdot,dmicro,"uvw-c-tpr-press-mdot-dmicro",0);
  }
  input.close();

#ifndef USE_VOF
  update_step(c, step, sflag);
#endif

  /* set iint (counter for TEC file) */
  int iint = int(time.current_time()/tint) + 1;
  boil::oout<<"iint= "<<iint<<"\n";

  const real rhol=liquid.rho()->value();
  const real rhov=vapor.rho()->value();

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

    /*-------------------+
    |  reset body force  |
    +-------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /*------------------+
    |  set heat source  |
    +------------------*/
    q=0.0;

    for_vk(tpr,k){
      if(approx(tpr.zn(k+1),0.0)){
        for_vij(tpr,i,j){
          if( fabs(tpr.xc(i))<0.5*h_width && fabs(tpr.yc(j))<0.5*h_width ){
            q[i][j][k]=qsrc*tpr.dV(i,j,k);  // [W/m3]*[m3] = [W]
          }
	}
      }
    }

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();    // calculate mdot in bulk
    pc.micro(&xyz); // calculate mdot in microlayer and update microlayer thickness
#ifndef USE_VOF
    update_step(c, step, sflag);  // need to update step for with IB
#endif
    ns.vol_phase_change(&f);  // outlet velocity is modified due to mdot

    /* bottom area */
    real area_l    = pc.get_hflux_area_l(Dir::ibody());
    real area_v    = pc.get_hflux_area_v(Dir::ibody());
    real area_micro= pc.get_hflux_area_micro(Dir::ibody());
    boil::oout<<"area= "<<time.current_time()<<" "<<area_l<<" "<<area_v
              <<" "<<area_micro<<"\n";

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* gravity force */
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k){
#ifdef USE_VOF
      real phil=c[i][j][k];
#else
      real phil=step[i][j][k];
#endif
      real phiv=1.0-phil;
      real deltmp=tpr[i][j][k]-tsat;
      real rhomix = (liquid.rho()->value() + liquid_drhodt*deltmp)*phil
                  + (vapor.rho()->value()  + vapor_drhodt*deltmp)*phiv;
      if(d.ibody().on(m,i,j,k))
        xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
    }

    /* surface tension */
#ifndef USE_VOF
    conc.tension(&xyz, mixed, step);
#else
    conc.tension(&xyz, mixed);
#endif
    //boil::plot->plot(xyz,c,mdot,"bodyforce-xyz-c-mdot",1);

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* turbulence model */
    //mod.smagorinsky( &ns, &mu_t, 0.173 );  // if this is commented out, then laminar flow

#if 1 
    /* increase viscosity in outlet region to stabilize calculation */
    const real z0=0.010;
    const real z1=0.015;
    for_avk(c,k){
      if(c.zc(k)>z0){
        real coef=std::min((c.zc(k)-z0)/(z1-z0),1.0);
        for_avij(c,i,j){
          mu_t[i][j][k]= coef * liquid.mu()->value() * 1000;
        }
      }
    }
#endif
    ns.discretize( &mu_t );
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step(); // cal convection term
    ns.grad(press);     // cal nabla.p
    ns.solve(ResRat(1e-14));  // cal diffusion term, and intermediate velocity

    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-6))) OMS(converged);
    p.exchange();
    ns.project(p);  // update velocity
    press += p;     // update pressure

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(d.ibody().on(i,j,k)){
        if(pmin>press[i][j][k]) pmin=press[i][j][k];
      }
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k){
      if(d.ibody().on(i,j,k)){
        press[i][j][k] -= pmin;
      } else {
        press[i][j][k] = 0.0;
      }
    }
    press.bnd_update();
    press.exchange_all();

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc.advance();
    conc.totalvol();

    /*---------------------------+
    |  replant seed or cut neck  |
    +---------------------------*/ 
    nucl.replant();

    /*-------------------------------+
    |  outlet region: delete bubbles |
    +-------------------------------*/
    /* delete bubble in outlet region: flash method */
    int kkm = -1;
    real clrmin_m=boil::exa;
    real clrmin_p=boil::exa;
    if (d.contains_z(z0)) {
      kkm=c.k(z0);
      for_vij(c,i,j) {
        if(clrmin_m>c[i][j][kkm-1]) clrmin_m=c[i][j][kkm-1];
      }
      for_vij(c,i,j) {
        if(clrmin_p>c[i][j][kkm]) clrmin_p=c[i][j][kkm];
      }
    }
    boil::cart.min_real(&clrmin_m);
    boil::cart.min_real(&clrmin_p);
    if ( (clrmin_m>0.5) && (clrmin_p>0.5) ){
      for_avk(c,k){
        if(c.zc(k)>z0){
          for_avij(c,i,j){
            c[i][j][k]= 1.0;
          }
        }
      }
    }
    /* update clr in cipcsl2 after seed, cutneck and outlet-region */
    c.bnd_update();
    c.exchange_all();
#ifndef USE_VOF
    conc.update_node(c);  // update color function at node, edge and face in accordance with c
    update_step(c, step, sflag);
#endif

    // temperature in outlet region
    for_avk(tpr,k){
      if(c.zc(k)>z0){
        for_avij(c,i,j){
          tpr[i][j][k]= tout;
        }
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();  // cal convection term
    enthFD.solve(ResRat(1e-16),"enthFD");  // cal diffusion term

    boil::oout<<"hflux:time,_total[W],_micro[W] "
              <<time.current_time()<<"  "
              <<pc.get_hflux(Dir::ibody())<<" "
              <<pc.get_hflux_micro(Dir::ibody())<<"\n";

    /*-------------+
    |  dt control  |
    +-------------*/
    /* minimum color function */
    conc.color_minmax();
    boil::oout<<"main:color_min,max= "<<time.current_time()<<" "
             <<conc.minval()<<" "<<conc.maxval()<<"\n";
    real clrmin=conc.minval();

    bool request_replant=false;
    for(int nsd=0; nsd<nucl.size(); nsd++){
      if( nucl.sites[nsd].req_replant() ) request_replant=true;
    }

    if ( clrmin < 0.5 || request_replant ) {
      /* interface is included */
      time.control_dt(ns.cfl_max(), cfl_with, dt);
    } else {
      /* interface is not included */
      time.control_dt(ns.cfl_max(), cfl_wo, 0.002);
    }

    real cflmax = ns.cfl_max();
    boil::oout<<"main:cflmax= "<<cflmax<<"\n";
    for(int nsd=0; nsd<nucl.size(); nsd++){
      if( cflmax<=cfl_with && time.dt()<=dt) {
        nucl.sites[nsd].set_allow_replant(true);
      } else {
        nucl.sites[nsd].set_allow_replant(false);
      }
    }

    boil::oout<<"main:request= "<<request_replant<<" "
              <<nucl.sites[0].allow_replant()<<" "
              <<cflmax<<" "<<time.dt()<<" "<<dt<<"\n";

    /*--------------+
    |  output data  |
    +--------------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      if(clrmin<0.5) {
        tpr.exchange_all();
        boil::plot->plot(uvw,c,tpr,press,mdot,dmicro
                       ,"uvw-c-tpr-press-mdot-dmicro",iint);
      }
      iint = int(time.current_time()/tint) + 1;
    }

    /* diameter */
    if(nucl.size()==1) {
      real dia=0.0;
      if ((area_v+area_micro)>0.0) {
        real zft=0.006;
        conc.front_minmax( Range<real>(-LX2,LX2) ,Range<real>(-LX2,LX2)
                          ,Range<real>(0, zft));
        dia = 0.5 * ( (conc.get_xmaxft() - conc.get_xminft())
                    + (conc.get_ymaxft() - conc.get_yminft()) );
      }
      boil::oout<<"Diameter= "<<time.current_time()<<" "<<dia
                <<" xcent "<<0.5*(conc.get_xminft()+conc.get_xmaxft())
                <<" ycent "<<0.5*(conc.get_yminft()+conc.get_ymaxft())<<"\n";
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % nint == 0 ) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
#ifndef USE_VOF
      conc .save("conc",  time.current_step());
#else
      c.save("c",  time.current_step());
#endif
      tpr  .save("tpr",   time.current_step());
      nucl .save("nucl",   time.current_step());
      dmicro.save("dmicro",time.current_step());
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
    if( boil::timer.current_min() > wmin
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
#ifndef USE_VOF
      conc .save("conc",  time.current_step());
#else
      c.save("c",  time.current_step());
#endif
      tpr  .save("tpr",   time.current_step());
      nucl .save("nucl",   time.current_step());
      dmicro.save("dmicro",time.current_step());
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
#ifndef USE_VOF
      conc .rm("conc", ts);
#else
      c.rm("c", ts);
#endif
      tpr  .rm("tpr", ts);
      nucl .rm("nucl", ts);
      exit(0); 
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}
