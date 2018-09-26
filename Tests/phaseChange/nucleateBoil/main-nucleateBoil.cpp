#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#include "update_step.cpp"

#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

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
const real h_width = 0.010;  // heater width

/* parameter for boundary or initial condition */
const real tsat   = 1.0E+02;
const real tout = tsat - 0.5;
const real twall = tsat + 9.0;
const real radius=0.0e-3;

/* constants */
const real gravity = 9.8;
const real pi = acos(-1.0);

/******************************************************************************/
main(int argc, char ** argv) {

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
  Grid1D gx0( Range<real>(-LX2,-LX1)
            , Range<real>(2.0*LX1/real(NX1),LX1/real(NX1))
            , NX2, Periodic::no() );
  Grid1D gx1( Range<real>(-LX1,LX1), 2*NX1, Periodic::no() );
  Grid1D gx2( Range<real>( LX1, LX2)
            , Range<real>(LX1/real(NX1),2.0*LX1/real(NX1))
            , NX2, Periodic::no() );
  Grid1D gt1 (gx0 , gx1, Periodic::no());
  Grid1D gx  (gt1 , gx2, Periodic::no());
  Grid1D gz0( Range<real>(LZ0, 0.0), NZ0, Periodic::no() );
  Grid1D gz1( Range<real>(0.0, LZ1), NZ1, Periodic::no() );
  Grid1D gz2( Range<real>(LZ1, LZ2)
            , Range<real>(LZ1/real(NZ1),4.0*LZ1/real(NZ1))
            , NZ2, Periodic::no() );
  Grid1D gztmp ( gz0, gz1, Periodic::no());
  Grid1D gz    ( gztmp, gz2, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Domain d(gx, gx, gz, & floor);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);           // vel
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d), cold(d); // color function
  Scalar tpr(d), q  (d), step(d), sflag(d); // temperature
  Scalar mdot(d);                  // phase-change rate
  Scalar dmicro(d);                // micro-layer film thickness
  Scalar mu_t(d);                  // eddy viscosity

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    //uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    //uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
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
  f = p.shape();
  mdot = p.shape();
  q = p.shape();
  mu_t = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  g = c.shape();
  step = c.shape();
  sflag = c.shape();
  dmicro = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), tout ) );

  /* heater power */
  real qsrc;
  real qflux=28.7*1000.0;     // [W/m2]
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
  vapor  .mu    (1.255e-5);
  vapor  .rho   (0.597);
  vapor  .cp    (2030*0.597);
  vapor  .lambda(0.025);
  liquid.mu    (0.28e-3);
  liquid.rho   (958.4);
  liquid.cp    (4215.9*958.4);
  liquid.lambda(0.679);
  sapphire.rho    (3980.0);
  sapphire.cp     (750.0*3980.0);
  sapphire.lambda (35.0);
  const real latent=2258.0*1e3;
  Matter mixed(liquid, vapor, & step);
  mixed.sigma(5.9e-2);
  const real liquid_drhodt=-0.7;   //[kg/m3K]
  const real vapor_drhodt=-0.0017; //[kg/m3K]

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 500000;
  const real tint = 1.0e-2;
  const real tint2 = 1.0e-2;
  const int  nint2= 1000*gLevel;
  const real dxmin = d.dxyz_min();
  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<"\n";
  const real dt  =10.0*pow(0.5*pow(dxmin,3.0)
                 / (2.0*3.1415*mixed.sigma()->value()),0.5);
  const real cfl_with = 0.05;
  const real cfl_wo   = 0.2;
  Times time(ndt, 0.002);
  time.print_time(false);
  time.set_coef_dec(0.75);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::ic2());
  Krylov * solver2 = new CG(d, Prec::di());

  /* color function */
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_itsharpen(10);
  conc.set_cangle(0.0);

  /* momentum equation */
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.convection_set(TimeScheme::forward_euler());

  /* pressure solver */
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  /* enthalpy equation */
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver2, & mixed ,tsat, & sapphire);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  /* nucleation site */
  real rseed=83.33333e-6;
  //real rseed=125.0e-6;
  //real rseed=187.5e-6;
  real area_neck=1.50e-7;
  if(gLevel>2) area_neck=1.50e-7;
  real zplant=-1.0; // when bottom of bubble reaches zplant, next seed is set
  real tseed =109.0;  // when temp of nucleation site reaches tseed, next ...
  const real dmicro_min=1.0e-8;
  Nucleation nucl( &c, &tpr, &q, &time, dmicro, &mixed, rseed,
                   dmicro_min, latent, conc.get_cangle());
  nucl.set_slope(1.0*4.46e-3);
  /* phase change */
  PhaseChange pc( mdot, tpr, q, c, g, f, step, uvw, time, &mixed , latent
                , tsat, &sapphire, &nucl);

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
    conc. load("conc", ts);
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
    for_vijk(c,i,j,k) 
      c[i][j][k] = 1.0;

    /* set seed */
    real zsite=rseed*cos(45.0/180.0*pi);
    nucl.add(Site( 0.000,  0.000, zsite, tseed, zplant));
#if 0
    /* plant seed from initial step */
    for(int ns=0; ns<nucl.size(); ns++){
      nucl.sites[ns].set_time_seed(0.0);
    }
    nucl.plant();
#endif
    c.bnd_update();
    c.exchange_all();

#if 0
    const real ztconst = 1.0*1.0e-3;
    for_vijk(c,i,j,k) {
      if(tpr.zc(k)<0.0) {
        tpr[i][j][k] = twall;
      } else if (tpr.zc(k)<ztconst) {
        tpr[i][j][k] = twall + (tout-twall)/ztconst * tpr.zc(k);
      } else {
        tpr[i][j][k] = tout;
      }
    }
#else
    tpr=tout;
#endif

    tpr.bnd_update();
    tpr.exchange_all();
    conc.init();
    boil::plot->plot(uvw,c,tpr,press,mdot,dmicro
                    ,"uvw-c-tpr-press-mdot-dmicro",0);
  }
  input.close();

  update_step(c, step, sflag);

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;
  int iint2 = int(time.current_time()/tint2) + 1;
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

    /*------------------+
    |  reset body force |
    +------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;
    q=0.0;

    for_vk(tpr,k){
      if(approx(tpr.zn(k+1),0.0)){
        for_vij(tpr,i,j){
          if( fabs(tpr.xc(i))<0.5*h_width && fabs(tpr.yc(j))<0.5*h_width ){
            q[i][j][k]=qsrc*tpr.dV(i,j,k);
          }
	}
      }
    }

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();
    pc.micro(&xyz);
    update_step(c, step, sflag);  // 0119 need to update step for with IB
    ns.vol_phase_change(&f);

    /* bottom area */
    real area_l=pc.get_hflux_area_l(Dir::ibody());
    real area_v=pc.get_hflux_area_v(Dir::ibody());
    boil::oout<<"area= "<<time.current_time()<<" "<<area_l<<" "<<area_v<<"\n";

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* gravity force */
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k){
      real phil=step[i][j][k];
      real phiv=1.0-phil;
      real deltmp=tpr[i][j][k]-tsat;
      real rhomix = (liquid.rho()->value() + liquid_drhodt*deltmp)*phil
                  + (vapor.rho()->value()  + vapor_drhodt*deltmp)*phiv;
      if(d.ibody().on(m,i,j,k))
        xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
    }

    /* surface tension */
    conc.tension(&xyz, mixed, step);
    //boil::plot->plot(xyz,c,mdot,"bodyforce-xyz-c-mdot",1);

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* essential for moving front */
    //mod.smagorinsky( &ns, &mu_t, 0.173 );
#if 1 
    /* increase viscosity in outlet region */
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
    ns.new_time_step();

    ns.grad(press);
    ns.solve(ResRat(1e-14));

    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-6))) OMS(converged);
    p.exchange();
    ns.project(p);
    press += p;

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
    cold = c;
    conc.advance();
    conc.totalvol();

    pc.modify_vel(uvw,c,cold);

    /*---------------------------+
    |  replant seed or cut neck  |
    +---------------------------*/ 
    nucl.replant();

    /*-------------------------------+
    |  outlet region: delete bubbles |
    +-------------------------------*/
#if 0
    /* delete gradually */
    //const real z0=0.011;
    //const real z1=0.015;
    for_avk(c,k){
      if(c.zc(k)>z0){
        real coef=std::min((c.zc(k)-z0)/(z1-z0),1.0);
        for_avij(c,i,j){
          c[i][j][k]= (1.0-coef)*c[i][j][k] + coef*1.0;
        }
      }
    }
#else
    /* flash */
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
#endif
    for_avk(tpr,k){
      if(c.zc(k)>z0){
        for_avij(c,i,j){
          tpr[i][j][k]= tout;
        }
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();

    /* update clr in cipcsl2 after seed, cutneck and outlet-region */
    c.bnd_update();
    c.exchange_all();
    conc.update_node(c);
    update_step(c, step, sflag);

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthFD");

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

#if 1
    if((time.current_time()) / (tint2) >= real(iint2) ) {
      iint2 = int(time.current_time() / tint2);
      std::cout.setf(std::ios_base::scientific);
      std::cout<< std::setprecision(16);
      /* temperature on wall */
      for_vk(c,k){
        if(approx(c.zn(k), 0.0, 1.0e-12)) {
          for_vi(c,i) {
            if(approx(c.xn(i), 0.0, 1.0e-12)) { 
              std::ofstream fout;
              std::stringstream ss;
              ss <<"twall_p"<<boil::cart.iam()<<"_"<<iint2<<".dat";
              std::string fname = ss.str();
              int len = fname.length();
              char * cfname = new char[len+1];
              memcpy(cfname, fname.c_str(), len+1);
              fout.open(cfname);
              for_vj(c,j) {
                real dw = 0.5 * c.dzc(k-1);
                real dflu;
                real lambdas = sapphire.lambda()->value();
                real lambdaf;
                real tflu;
                if( dmicro[i][j][k]>boil::mega ){
                  dflu = 0.5 * c.dzc(k);
                  lambdaf=liquid.lambda()->value();
                  tflu = tpr[i][j][k];
                } else if( dmicro[i][j][k]<dmicro_min+boil::pico ){
                  dflu = 0.5 * c.dzc(k);
                  lambdaf=vapor.lambda()->value();
                  tflu = tpr[i][j][k];
                } else {
                  dflu=dmicro[i][j][k];
                  lambdaf=liquid.lambda()->value();
                  tflu=tsat;
                }
                  
                real tw = ( dflu * lambdas * tpr[i][j][k-1]
                          + dw * lambdaf * tflu )
                        / ( dflu * lambdas + dw * lambdaf);

                fout<<c.xc(i)<<" "<<c.yc(j)<<"  "<<tw<<" "<<tpr[i][j][k-1]<<" "
                    <<tpr[i][j][k]<<" "<<dmicro[i][j][k]<<" "
                    <<c[i][j][k]<<" "<<i<<" "<<j<<" "<<k<<"\n";
              }
              fout.close();
            }
          }
        }
      }
      std::cout<< std::setprecision(8);
#endif
      iint2 = int(time.current_time()/tint2) + 1;
    }

    /* diameter */
    if(nucl.size()==1) {
      real dia=0.0;
      if (area_v>0.0) {
        real zft=0.006;
        conc.front_minmax( Range<real>(-LX2,LX2) ,Range<real>(-LX2,LX2)
                          ,Range<real>(0, zft));
        dia = 0.5 * ( (conc.get_xmaxft() - conc.get_xminft())
                    + (conc.get_ymaxft() - conc.get_yminft()) );
      }
      boil::oout<<"Diameter= "<<time.current_time()<<" "<<dia<<"\n";
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if((time.current_step()) % (nint2)==0 ) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
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
    if( boil::timer.current_min() > (wmin-30.0)
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
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
      conc .rm("conc", ts);
      tpr  .rm("tpr", ts);
      nucl .rm("nucl", ts);
      exit(0); 
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}

/*-----------------------------------------------------------------------------+
 '$Id: main-nucleateBoil.cpp,v 1.1 2018/04/30 08:48:02 sato Exp $'/
+-----------------------------------------------------------------------------*/
