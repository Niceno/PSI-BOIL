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
const int gLevel = 2;  //grid level=2,3,4
const real LX2 =  0.15 *0.001;
const real LX1 =  0.1*0.001;
const int NX1  = 25*gLevel;
const int NX2  =  7*gLevel;
const real dxuni = LX1/real(NX1);
const real LZ0 = -8.0e-6;
const real LZ1  = 0.1*0.001;
const real LZ2  = 0.15*0.001;
const int NZ0  =  4;  // solid
const int NZ1  = 25*gLevel;
const int NZ2  =  7*gLevel -4;

/* parameter for boundary or initial condition */
//const real tsat = 257.44+273.15;
const real tsat = 530.18;
const real tout = tsat;
const real twall = tsat + 1.822;

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
  Grid1D gx1( Range<real>(0.0,LX1), NX1, Periodic::no() );
  Grid1D gx2( Range<real>(LX1,LX2),
              Range<real>(dxuni,3.0*dxuni),
              NX2, Periodic::no() );
  Grid1D gx ( gx1, gx2, Periodic::no());
  Grid1D gz0( Range<real>(LZ0, 0.0), NZ0, Periodic::no() );
  Grid1D gz1( Range<real>(0.0, LZ1), NZ1, Periodic::no() );
  Grid1D gz2( Range<real>(LZ1, LZ2),
              Range<real>(1.2*dxuni, 3.0*dxuni),
              NZ2, Periodic::no() );
  Grid1D gztmp ( gz0, gz1, Periodic::no());
  Grid1D gz ( gztmp, gz2, Periodic::no());

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
  Scalar c  (d), g  (d), kappa(d); // concentration
  Scalar tpr(d), q  (d);           // temperature
  Scalar step(d), sflag(d);        // step function
  Scalar mdot(d);                  // phase-change rate
  Scalar dmicro(d);                // micro-layer film thickness
  Scalar mu_t(d);                  // eddy viscosity

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
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

  c.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  g = c.shape();
  step = c.shape();
  sflag = c.shape();
  dmicro = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(),531.091 ) );

  /* heater power */
  real qflux=7.29*1000.0;     // [W/m2]
  real qsrc=qflux/fabs(LZ0);  // [W/m3]
  boil::oout<<"#qsrc= "<<qsrc<<"\n";

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d), sapphire(d);
  vapor  .mu    (1.78e-5);
  vapor  .rho   (22.539);
  vapor  .cp    (4214.5*22.539);  // J/m3
  vapor  .lambda(0.053181);
  liquid.mu    (1.03e-4);
  liquid.rho   (788.25);
  liquid.cp    (4949.3*788.25);   // J/m3
  liquid.lambda(0.61293);
  sapphire.rho    (8908.0);//changing to nickel ??
  sapphire.cp     (444.0*8908.0);
  sapphire.lambda (90.9);
  const real latent=1677.9*1000.0;
  Matter mixed(liquid, vapor, &step);
  mixed.sigma(0.024);
  const real liquid_drhodt=-1.6;   //[kg/m3K]
  const real vapor_drhodt=-0.101; //[kg/m3K]

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 700000;
  const real tint = 1.0e-4;
  const int  nint2= 1000;
  const real dxmin = d.dxyz_min();
  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<"\n";
  const real dt  =10.0*pow(vapor.rho()->value()*pow(dxmin,3.0)
                 / (2.0*3.1415*mixed.sigma()->value()),0.5);
  const real cfl_with = 0.05;
  const real cfl_wo   = 0.2;
  Times time(ndt, dt);
  //Times time(ndt, 0.002);
  time.print_time(false);
  time.set_coef_dec(0.75);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::ic2());
  Krylov * solver2 = new CG(d, Prec::di());

  /* color function */
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(10);
  conc.set_cangle(0.0);
  conc.set_globalSharpen();

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
  real rseed = 4.0e-6;
  real zplant= 0.05*0.001; // when bottom of bubble reaches zplant, next seed is set
  real tseed = twall-0.001;  // when temp of nucleation site reaches tseed, next ...
  //Nucleation nucl( &c, &tpr, &time, dmicro, rseed
  //               , mixed.sigma()->value(), 0.0e-10, conc.get_cangle());
  Nucleation nucl( &c, &tpr, &q, &time, dmicro, &mixed, rseed
                 , 0.0e-10, latent, conc.get_cangle());
  nucl.set_slope(1.0*4.46e-3);
  nucl.set_seed_period(0.001);
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
    real zsite=rseed*cos(0.0/180.0*pi);
    nucl.add(Site( 0.000,  0.000, zsite, tseed, zplant));
#if 1
    /* plant seed from initial step */
    nucl.plant();
    for(int ns=0; ns<nucl.size(); ns++){
      nucl.sites[ns].set_time_seed(0.0);
    }
#endif
    c.bnd_update();
    c.exchange_all();

    const real ztconst = 0.324*1.0e-3;
    for_vijk(c,i,j,k) {
      if(tpr.zc(k)<0.0) {
        tpr[i][j][k] = twall;
      } else {
        tpr[i][j][k] = twall + ((tout+0.911)-twall)/ztconst * tpr.zc(k);
      } 
    }
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

    /*----------------------------------------------+
    |  reset body force & source term for enthalpy  |
    +----------------------------------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    q=0.0; 
    for_vk(tpr,k){
      if((tpr.zc(k)<0.0)){
        for_vij(tpr,i,j){
          q[i][j][k]=qsrc*tpr.dV(i,j,k);
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
    const real z0=0.25*0.001;//JSM:changed from0.3
    const real z1=0.324*0.001;//JSM:changed from 0.4

    for_avk(c,k){
      if(c.zc(k)>z0){
        real coef=std::min((c.zc(k)-z0)/(z1-z0),1.0);
        for_avij(c,i,j){
          mu_t[i][j][k]= coef * liquid.mu()->value() * 10;
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
    //pr.solve(ResRat(1e-6));
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
    conc.advance();

    /*---------------------------+
    |  replant seed or cut neck  |
    +---------------------------*/ 
    nucl.replant();

    /*-------------------------------+
    |  outlet region: delete bubbles |
    +-------------------------------*/
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
    //std::cout<<"main:clrmin_m= "<<clrmin_m<<" "<<clrmin_p<<" "<<kkm<<"\n";
    if ( (clrmin_m>0.5) && (clrmin_p>0.5) ){
      for_avk(c,k){
        if(c.zc(k)>z0){
          for_avij(c,i,j){
            c[i][j][k]= 1.0;
          }
        }
      }
    }

//JSM:commented out because TBL is greater than domain size
//added if condition
#if 0
    for_avk(tpr,k){
      if(c.zc(k)>z0){
        for_avij(c,i,j){
          tpr[i][j][k]= tout;
        }
      }
    }
#endif
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

#if 0
    /* enforce twall in solid domain */
    for_vk(tpr,k){
      if(tpr.zc(k)<0.0){
        for_vij(tpr,i,j){
          tpr[i][j][k]=twall;
        }
      }
    }
#endif

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
      tpr.exchange_all();
      boil::plot->plot(uvw,c,tpr,press,mdot,dmicro
                     ,"uvw-c-tpr-press-mdot-dmicro",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /* diameter */
    if(nucl.size()==1) {
      real dia=0.0;
      real height=0.0;
      if (area_v>0.0) {
        real zft=0.006;
        conc.front_minmax( Range<real>(-LX2,LX2) ,Range<real>(-LX2,LX2)
                          ,Range<real>(0, zft));
        dia = 0.5 * ( (conc.get_xmaxft() - conc.get_xminft())
                    + (conc.get_ymaxft() - conc.get_yminft()) );
        height = conc.get_zmaxft();
      }
      boil::oout<<"Diameter= "<<time.current_time()<<" "<<dia<<" "
                <<height<<"\n";
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if((time.current_step()) % (nint2)==0 ){
     //|| time.current_step() == 18273 ) {
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
        //output.open("time.txt", std::ios::out);
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
 '$Id: main-sakashita.cpp,v 1.1 2018/04/30 08:52:20 sato Exp $'/
+-----------------------------------------------------------------------------*/
