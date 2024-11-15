#include "Include/psi-boil.h"
#include <fstream>
#include <iostream>
#include <fenv.h>

#define VARIABLE

/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc<2){
    boil::oout<<"Two command line arguments required!"<<"\n";
    boil::oout<<"./Boil wmin"<<"\n";
    exit(0);
  }

/******************************************************************************/
/* ------------ input from command line */
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  const int gLevel = 1;
  boil::oout<<"glevel= "<<gLevel<<"\n";

/******************************************************************************/
/* ------------ boundary or initial conditions */
  const real velin = -0.1;   // velocity inlet [m/s] for both phases
  const real tref  = 373.15; // reference temperature [K], =Tsat(eps=0)
                             // Don't change this unless system prssure
                             // or working fluid is changed!
  const real tin   = 373.15; // inlet temperature for both phases
  const real tout  = 363.15; // outlet temperature
  const real tinit = 373.15; // initial temperature
  const real eps0  = 0.3000; // eps: volume fraction of air in gas phase [-]
                             // used for the inlet and initial condition
  const real gravity = boil::g;  // [m/s2]

/******************************************************************************/
/* ------------ numerical simulation settings */
  /* total number of steps */
  const int ndt = 35000*gLevel; /* inconsequential */

  /* total time */
  const real tend = 0.5;

  /* steps per backup */
  const int n_per_backup= 5000;

  /* plotting each t_per_plot seconds */
  const real t_per_plot = 2.0e-4;

  /* dt settings */
  const real surftens_dt_coef = 0.126;
  const real initdtcoef = 1e-1;

  /* cfl limit */
  const real cfl_limit = 0.1;//0.15;//0.05;

/* ------------ optional simulation settings */

  /* multigrid */
  const bool multigrid_stop_if_diverging = true;
  //const bool multigrid_stop_if_diverging = false;

  const int multigrid_min_cycles = 1;
  const int multigrid_max_cycles = 20;

  ResTol rt = ResTol(-1.);
  ResRat rr = ResRat(1e-5);

  MaxIter multigrid_mm_smooth1 = MaxIter(35);
  MaxIter multigrid_mm_smooth2 = MaxIter(40);
  MaxIter multigrid_mm_solve = MaxIter(110);
  MaxIter multigrid_mm_stale1 = MaxIter(-1);
  MaxIter multigrid_mm_stale2 = MaxIter(-1);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm_smooth1,multigrid_mm_smooth2,multigrid_mm_solve};
  std::array<MaxIter,3> multigrid_mstale = {multigrid_mm_stale1,multigrid_mm_stale1,multigrid_mm_stale2};

  /* vof */
  const CurvMethod curv_method = CurvMethod::HF();

  const bool use_fs_interp = false;
  const bool store_pressure_extrap = true;
  const int niter_pressure_extrap = 8000;

  /* phase change - 4 version */
  const AccuracyOrder pc_accord = AccuracyOrder::FourthUpwind();
  const bool discard_points_near_interface = false;
  const bool use_unconditional_extrapolation = false;

  /* under-relaxation */
  const real ur = 0.01;

  /* limit vfrac for ctp */
  const real limitvf = 0.001;

  /* moving-frame */
  const real ur_vel = 0.1;
  real wvel = 0.0;

/******************************************************************************/
/* ------------ material properties */
  const real Mv = 18.0e-3;  // molecular mass of water [kg/mol] used for diffcoef and clap

  /* gas modelled as dry air */
  const real mug = 1.962e-5;
  const real rhog = 1.086;
  const real cpg = 1.0063e3*rhog;
  const real lambdag = 2.816e-2;
  const real diffcoef = 0.26e-4; // diffusion coefficient of air-steam [m2/s] at 15 deg

  // water at 20 deg
  const real mul = 1.00e-3;
  const real rhol = 998.2;
  const real cpl = 4181.8*rhol;
  const real lambdal = 0.5984;

  real sig = 0.07274;
  const real latent=2.0e6; 

  const real betal = 0.;
  const real betag = 0.;//1./tsat0; /* ideal gas approximation

/******************************************************************************/
/* ------------ domain dimensions */
  const real thick_film = 2.0e-4;   // liquid film thickness [m]
  const real LX1 = 1.5*thick_film;  // fine grid will be made in 0 <= x <= LX1
  const int  NX1 = 16*gLevel;       // no. cells for fine grid region
  const int  NX_film = 16*gLevel/1.5;
  const real DX  = LX1/real(NX1);
  boil::oout<<"NX_film,DX= "<<NX_film<<" "<<DX<<"\n";

  const int   NX2 = 16*gLevel;
  const real  LX2 = 8.0 * LX1;
  const int   NX = NX1 + NX2;

  const real DZ = 4.0*DX;
  const int  NZ = 32 * gLevel;
  const real LZ = real(NZ)*DZ;

  const int  NY = 8 *gLevel;
  const real LY = real(NY)*DZ;

/******************************************************************************/
/* ------------- setup finished */
/******************************************************************************/
/* below this line, NO VALUES AND SETTINGS CAN BE ENTERED! */
/******************************************************************************/

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx0( Range<real>(0.0, LX1), NX1, Periodic::no() );
  Grid1D gx1( Range<real>(LX1, LX2), Range<real>(1.2*DX,16.0*DX), 
              NX2, Periodic::no() );
  //Grid1D gx (gx0, gx1, Periodic::no());
  Grid1D gx (gx0, gx1, Periodic::no(), BndGrid::wall(), BndGrid::symmetry());

  Grid1D gy ( Range<real>(-0.5*LY,0.5*LY) ,NY ,Periodic::yes() );
  Grid1D gz ( Range<real>(0.0,LZ) ,NZ ,Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);
  //boil::plot->plot(d,"domain");
  const real dxmin = d.dxyz_min();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);     // uvw: velocity [m/s], xyz: source term for NS (body force) [N]
  Vector uvw_1(d), uvw_2(d); // phasic velocity [m/s]
  Scalar p(d), f(d), press(d,"press"); // p: pressure difference in time step [Pa]
                                       // f: source term of pressure Poisson equation [m3/s2]
				       // press: pressure [Pa]
  Scalar c  (d,"vfl"), g  (d), kappa(d); // c: color function (= volume fraction of liquid) [-]
                                         // g: source term of the equation for color [1/s]
                                         // kappa: curvature [1/m]
  Scalar tpr(d,"T"), q  (d);           // tpr: temperature [deg.], q: heat source [W]
  Scalar mdot(d,"mdot"), mflx(d);      // mdot: phase-change rate [kg/m3s]
  Scalar eps(d,"eps"), mdot_eps(d);    // eps: volume fraction of air in gas phase [-]

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd(Range<int>(1,NX_film), Range<int>(1,NY), Dir::kmax(),
              BndType::inlet(),0.0, 0.0, velin ) );
    uvw.bc(m).add( BndCnd(Range<int>(NX_film+1,NX), Range<int>(1,NY), Dir::kmax(),
              BndType::inlet(),0.0, 0.0, velin ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw_1(m)=uvw(m).shape();
    uvw_2(m)=uvw(m).shape();
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  /* copy b.c. from p */
  press = p.shape();
  f = p.shape();
  mdot = p.shape();
  mflx = p.shape();
  q = p.shape();
  g = p.shape();
  kappa = p.shape();
  mdot_eps = p.shape();

  c.bc().add( BndCnd(Range<int>(1,NX_film), Range<int>(1,NY), Dir::kmax(),
              BndType::dirichlet(),1.0 ) );
  c.bc().add( BndCnd(Range<int>(NX_film+1,NX), Range<int>(1,NY), Dir::kmax(),
              BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  tpr.bc().add( BndCnd( Dir::imin(), BndType::neumann()));
  tpr.bc().add( BndCnd( Dir::imax(), BndType::symmetry()));
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), tin ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic()));
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic()));

  eps.bc().add( BndCnd( Dir::imin(), BndType::neumann()));
  eps.bc().add( BndCnd( Dir::imax(), BndType::symmetry()));
  eps.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), eps0 ) );
  eps.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), eps0 ) );
  eps.bc().add( BndCnd( Dir::jmin(), BndType::periodic()));
  eps.bc().add( BndCnd( Dir::jmax(), BndType::periodic()));

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter gas(d), liquid(d);
  gas  .mu    (mug);
  gas  .rho   (rhog);
  gas  .cp    (cpg);  // J/m3
  gas  .lambda(lambdag);
  gas.mmass    (Mv);
  gas.gamma    (diffcoef*rhog);
  gas.beta     (betag);
  liquid.mu    (mul);
  liquid.rho   (rhol);
  liquid.cp    (cpl);   // J/m3
  liquid.lambda(lambdal);
  liquid.beta  (betal);

  Matter mixed(liquid, gas, & c);
  mixed.sigma(sig);
  mixed.latent(latent);

  Matter * zero = &gas;
  Matter * one = &liquid;

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const real dt = surftens_dt_coef
                * pow(std::max(gas.rho()->value(),liquid.rho()->value())*pow(dxmin,3.0)
                / mixed.sigma()->value(),0.5);

  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<" "<<dt<<"\n";
  boil::oout<<"main:nparams= "<<t_per_plot<<" "<<n_per_backup<<"\n";
  Times time(ndt, dt);
  time.set_coef_dec(0.75);
  time.set_dt(dt*initdtcoef);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solverCGic = new CG(d, Prec::ic2());
  Krylov * solverCGdi = new CG(d, Prec::di());
  Krylov * solverBICG = new BiCGS(d, Prec::di());

  Krylov * solver = solverCGic;

  /*-------------------+
  |  define equations  |
  +-------------------*/
  /* momentum equation */
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.convection_set(TimeScheme::forward_euler());
  ns.diffusion_set(TimeScheme::backward_euler());

  /* pressure solver */
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(multigrid_stop_if_diverging);
  multigrid.min_cycles(multigrid_min_cycles);
  multigrid.max_cycles(multigrid_max_cycles);

  /* color function */
  Vector & uvwvof = uvw_1;
  VOF conc(c, g, kappa, uvwvof, time, solver);
  conc.set_curv_method(curv_method);
  conc.set_use_interp(use_fs_interp);
  conc.set_pressure_extrapolation_parameters(store_pressure_extrap,niter_pressure_extrap);
  //conc.set_advection_method(AdvectionMethod::ReconstructedSplit());

  /* enthalpy equation */
#ifndef VARIABLE
  TIF tsat(tref);
#else
  // Antoine model. saturation temperature in a cell, which does not include liquid-vapor 
  // interface is set to tref, but it is not used.
  Antoine tsat(tref,conc.topo,eps,8.07131,1730.63,233.426);
  tsat.set_ur(ur);
  real epsinf = tsat.epsilon(tout);
  real epstest = tsat.epsilon(tref);
  real volrat = eps0/epsinf;
  real radrat = pow(volrat,1./3.);
  boil::oout<<"Tint(eps="<<eps0<<")= "<<tsat.temperature(eps0)<<" [K] "
            <<tsat.temperature(eps0)-273.15<<" [C]\n";
  boil::oout<<"eps(Tsat="<<tout<<")= "<<epsinf<<"\n";
#endif

  CommonHeatTransfer cht(tpr,conc.topo,tsat,&mixed);

  Vector & uvwenth1 = uvw_1;
  Vector & uvwenth2 = uvw_2;
  EnthalpyFD enthFD (tpr, q, uvw, uvwenth1, uvwenth2, time, solver  , &mixed, cht);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  /* phase change */
  PhaseChange4 pc(mdot, mflx, q, g , f , uvw, cht, time, &mixed);
  pc.set_accuracy_order(pc_accord);
  pc.set_discard_points_near_interface(discard_points_near_interface);
  pc.set_unconditional_extrapolation(use_unconditional_extrapolation);

  /* ng transport */
  ConcentrationTP ngtransp(eps,mdot_eps,uvw,
                           conc.flow(),conc.heaviside(),conc.topo,
                           time, solverCGdi, &gas,limitvf);

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  boil::test_irun();
  boil::set_irun(1);

  int ts;
  /* load variables */
  std::vector<Scalar*> load_scalars = { &press, &c, &tpr, &eps };
  std::vector<std::string> load_scalar_names = { "press", "c", "tpr", "eps" };
#ifdef VARIABLE
  load_scalars.push_back(&(tsat.tif));
  load_scalar_names.push_back("tif");
#endif
  std::vector<Vector*> load_vectors = { &uvw };
  std::vector<std::string> load_vector_names = { "uvw" };

  if(boil::load_backup("time.txt",ts,time,
                       load_scalars, load_scalar_names,
                       load_vectors, load_vector_names)) {
    conc.init();
#ifdef VARIABLE
    tsat.update_tifold();
#endif

  } else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    /*--------------------+
    |  initial condition  |
    +--------------------*/
    Comp m = Comp::w();
    for_vmijk(uvw,m,i,j,k){
      uvw[m][i][j][k] = velin;
    }
    uvw.exchange();

    for_vijk(c,i,j,k) {
      if(c.xc(i)<thick_film){
        c[i][j][k] = 1.0;
      }
    }
    c.bnd_update();
    c.exchange_all();
    conc.init();
    conc.totalvol();

    tpr = tinit;
    tpr.bnd_update();
    tpr.exchange_all();

    eps = eps0;
    eps.bnd_update();
    eps.exchange_all();

#ifdef VARIABLE
    tsat.init();
#endif

    boil::plot->plot(uvw,c,tpr,eps,press,mdot,"uvw-c-tpr-eps-press-mdot",0); 
  }

  /* set iint */
  int iint = int(time.current_time()/t_per_plot) + 1;
  boil::oout<<"iint= "<<iint<<"\n";
  
  /*----------------+
  |  topology init  |
  +----------------*/
  conc.new_time_step();

  //real zpos(0.0);
  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /* temperature field */
    tsat.tint_field();

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();                // calculate mdot from temperature field
    // outlet region for mdot
    for_vijk(mdot,i,j,k){
      if (mdot.zc(k)<0.05*LZ) {
        mdot[i][j][k]=0.0;
      }
    }
    mdot.bnd_update();
    mdot.exchange_all();

    ngtransp.mdot_cutoff(mdot); // set mdot=0 for the cells of vfv < limitvf
    ns.vol_phase_change(&f);    // calculate volume change in whole domain

    /* reset body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* surface tension */
    conc.tension(&xyz, mixed,conc.color());
    // outlet
    for_m(m){
      for_vmijk(xyz,m,i,j,k){
        if(xyz.zc(m,k)<0.2*LZ){
          xyz[m][i][j][k]=0.0;
	}
      }
    }

    /* gravity force */
    Comp m = Comp::w();
    for_avmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
    }

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* essential for moving front */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();

    ns.grad(press);
    ns.solve(ResRat(1e-14));

    p = 0.;
    if(multigrid.cycle(Cycle::Z(),Cycle::F(),rt,rr,multigrid_mi,multigrid_mstale)) OMS(converged);
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

    ns.cfl_max();

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc.new_time_step();
    conc.advance_with_extrapolation(true,ResTol(1e-7),uvw,f,
                                    one,&uvw_1,zero,&uvw_2);

    // outlet region
    for_vijk(c,i,j,k){
      if (c.zc(k)<0.05*LZ) {
        c[i][j][k]=0.0;
      }
    }
    c.bnd_update();
    c.exchange_all();

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthFD");

    /*-------------------------+
    |  solve species equation  |
    +-------------------------*/
    ngtransp.discretize();
    ngtransp.new_time_step();
    ngtransp.solve(ResRat(1e-14),"Concentration");
    ngtransp.extrapolate();
    eps.bnd_update();
    eps.exchange_all();

    /*-------------+
    |  dt control  |
    +-------------*/
    /* minimum color function */
    conc.color_minmax();

    /* front */
    conc.front_minmax();

    time.control_dt(ns.cfl_max(),cfl_limit,dt);

    /*--------------+
    |  output data  |
    +--------------*/
    if((time.current_time()) / (t_per_plot) >= real(iint) || time.current_step()==1 ) {
      boil::plot->plot(uvw,c,tpr,eps,press,mdot,"uvw-c-tpr-eps-press-mdot",iint);
      iint++;
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % n_per_backup == 0) {
      boil::save_backup(time.current_step(), 0, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
    }

    if( boil::timer.current_min() > wmin
      || time.current_step()==time.total_steps()
      || time.current_time()>tend) {
      boil::oout<<"main:current_min= "<<boil::timer.current_min() <<" "<< wmin<<"\n";
      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names);

      if(boil::timer.current_min() > wmin) {
        boil::set_irun(0);
      }
      break;
    }

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
