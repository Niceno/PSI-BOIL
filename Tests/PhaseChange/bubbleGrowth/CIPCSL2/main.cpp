#include "Include/psi-boil.h"
#include <fstream>
#include <iostream>
#include <fenv.h>
#define _GNU_SOURCE 1
#if 1
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif

#define SYM
//#define FS_TEST
//#define USE_MG
#define CASE 4
/* case 1 = axisymmetric VOF
        2 = cartesian 2D VOF
        3 = cartesian 3D VOF
        4 = cartesian 3D CIPCSL2
*/

#if CASE == 4
  void set_stepfunc(Scalar & sca, Scalar & scb, Scalar & scc);
#endif

/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc<3){
    boil::oout<<"Two command line arguments required!"<<"\n";
    boil::oout<<"./Boil wmin glevel"<<"\n";

    exit(0);
  }

/******************************************************************************/
/* ------------ input from command line */
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  const int gLevel = atoi(argv[2]); /* domain dimensions */ 
  boil::oout<<"glevel= "<<gLevel<<"\n";

/******************************************************************************/
/* ------------ boundary or initial conditions */
  const real tsat0 = 373.15;
  const real tout = tsat0+1.25;

/******************************************************************************/
/* ------------ numerical simulation settings */

  /* total time */
  //const real tend = 0.003;
  const real tend = 0.0012*1.25; /* 12247 for level 2 */

  /* number of backup points */
  const int n_bck = 2;

  /* number of plot points */
  const int n_plot = 25;

  /* in-time step iteration */
  const int mSimple = 1;  // Lubomir's original = 2

  /* dt settings */
  const real surftens_dt_coef = 10.;

  /* cfl limit */
  const real cfl_limit = 0.1;

/* ------------ optional simulation settings */

  /* multigrid */
  const bool multigrid_stop_if_diverging = true;
  //const bool multigrid_stop_if_diverging = false;

  const int multigrid_min_cycles = 1;
  const int multigrid_max_cycles = 20;

  /* vof */
  const CurvMethod curv_method = CurvMethod::HF();

  const bool use_fs_interp = false;
  const bool store_pressure_extrap = true;
  const int niter_pressure_extrap = 1000;

/******************************************************************************/
/* ------------ material properties */
  const real muv = 1.255e-5;
  const real rhov = 0.597;
  const real cpv = 2030*rhov;
  const real lambdav = 0.025;

  const real mul = 0.28e-3;
  const real rhol = 958.4;
  const real cpl = 4.2159e3*rhol;
  const real lambdal = 679e-3;

  const real sig = 59e-3;
  const real latent = 2.258e6;

/******************************************************************************/
/* ------------ domain dimensions */
  const int NX = 24*gLevel;
#ifdef SYM
  const int NZ = 24*gLevel;
#else
  const int NZ = 48*gLevel;
#endif

  const real LX = 187.5e-6;
  const real LZ = LX;

  const real DX = LX/real(NX);
  const real radius=50.0e-6;

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
  Grid1D gx( Range<real>(0.0,LX), NX, Periodic::no(), BndGrid::symmetry(), BndGrid::wall() );
#ifndef SYM
  Grid1D gz( Range<real>(-LZ,LZ), NZ, Periodic::no() );
#else
  Grid1D gz( Range<real>(0.0,LZ), NZ, Periodic::no(), BndGrid::symmetry(), BndGrid::wall() );
#endif

  /*---------+
  |  domain  |
  +---------*/
#if CASE == 1
  Axisymmetric d(gx, gz, DX);
#elif CASE == 2
  Grid1D gy(DX);
  Domain d(gz,gy,gz);
#else
  Domain d(gz,gz,gz);
#endif
  const real dxmin = d.dxyz_min();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), uvw_old(d),xyz(d);// vel
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d); // concentration
  Scalar tpr(d), q  (d);           // temperature
  Scalar mdot(d), mflx(d);         // phase-change rate
  Scalar tprold(d), tprap1(d), tprap2(d);
#if CASE == 4
  Scalar step(d), sflag(d);
#else
  Vector uvw_1(d), uvw_2(d);       // phasic vel
#endif

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
#if CASE == 1
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
#else
  #ifdef SYM
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
  #else
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
  #endif
#endif
#ifdef SYM
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
#else
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
#endif
#if CASE < 3
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
#else 
  #ifdef SYM
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::outlet() ) );
  #else
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::outlet() ) );
  #endif
#endif
#if CASE != 4
    uvw_1(m)=uvw(m).shape();
    uvw_2(m)=uvw(m).shape();
#endif
    uvw_old(m)=uvw(m).shape();
  }

#if CASE == 1
  c.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
#else
  #ifdef SYM
  c.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  #else
  c.bc().add( BndCnd( Dir::imin(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  #endif
#endif
#ifdef SYM
  c.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
#else
  c.bc().add( BndCnd( Dir::kmin(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
#endif
#if CASE < 3
  c.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
#else
  #ifdef SYM
  c.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::outlet() ) );
  #else
  c.bc().add( BndCnd( Dir::jmin(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::outlet() ) );
  #endif
#endif

  press=c.shape();
  p=c.shape();
  kappa=c.shape();
  f = p.shape();
  mdot = p.shape();
  mflx = p.shape();
  q = p.shape();
  g = c.shape();
#if CASE == 4
  step = c.shape();
  sflag = c.shape();
#endif

#if CASE == 1
  tpr.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
#else
  #ifdef SYM
  tpr.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
  #else
  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
  #endif
#endif
#ifdef SYM
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), tout ) );
#else
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), tout ) );
#endif
#if CASE < 3
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) ); 
#else
  #ifdef SYM
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), tout ) );
  #else
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), tout ) );
  #endif
#endif
  tprold = tpr.shape();
  tprap1 = tpr.shape();
  tprap2 = tpr.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d);
  vapor  .mu    (muv);
  vapor  .rho   (rhov);
  vapor  .cp    (cpv);  // J/m3
  vapor  .lambda(lambdav);
  liquid.mu    (mul);
  liquid.rho   (rhol);
  liquid.cp    (cpl);   // J/m3
  liquid.lambda(lambdal);

#if CASE == 4
  Matter mixed(liquid, vapor, & step);
#else
  Matter mixed(liquid, vapor, & c);
#endif
  mixed.sigma(sig);
  mixed.latent(latent);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const real dt = surftens_dt_coef*pow(vapor.rho()->value()*pow(dxmin,3.0)
                / (2.0*boil::pi*mixed.sigma()->value()),0.5);
  const int ndt = tend/dt;
  const int n_per_backup = ndt/n_bck;
  const int n_per_plot = ndt/n_plot;

  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<" "<<dt<<"\n";
  boil::oout<<"main:nparams= "<<n_per_plot<<" "<<n_per_backup<<"\n";
  Times time(ndt, dt);
  time.set_coef_dec(0.75);
  time.set_dt(dt);

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
  if(mSimple>1)
    ns.convection_set(TimeScheme::backward_euler()); //ns.convection is mandatory
  else
    ns.convection_set(TimeScheme::forward_euler());
  ns.diffusion_set(TimeScheme::backward_euler());

  /* pressure solver */
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(multigrid_stop_if_diverging);
  multigrid.min_cycles(multigrid_min_cycles);
  multigrid.max_cycles(multigrid_max_cycles);

  /* color function */
#if CASE != 4
  Vector & uvwvof = uvw_1;
#endif
#if CASE > 3
  boil::oout<<"main:use CIPCSL2\n";
  CIPCSL2 conc(c, g, kappa, uvw, time, solver);
#elif CASE > 1
  boil::oout<<"main:use VOF\n";
  VOF conc(c, g, kappa, uvwvof, time, solver);
#else
  boil::oout<<"main:use VOFaxisym\n";
  VOFaxisym conc(c, g, kappa, uvwvof, time, solver);
#endif
#if CASE < 4
  conc.set_curv_method(curv_method);
  conc.set_use_interp(use_fs_interp);
  conc.set_pressure_extrapolation_parameters(store_pressure_extrap,niter_pressure_extrap);
#endif

  /* enthalpy equation */
  TIF tsat(tsat0);
  CommonHeatTransfer cht(tpr,conc.topo,tsat,&mixed);

#if CASE != 4
  Vector & uvwenth1 = uvw_1;
  Vector & uvwenth2 = uvw_2;
#endif
#if CASE > 3
  EnthalpyFD enthFD (tpr, q, uvw, time, solver, &mixed, cht);
#elif CASE > 1
  EnthalpyFD enthFD (tpr, q, uvw, uvwenth1, uvwenth2, time, solver, &mixed, cht);
#else
  EnthalpyFDaxisym enthFD(tpr, q, uvw, uvwenth1, uvwenth2, time, solver  , &mixed,
                          cht);
#endif
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  enthFD.set_flux_accuracy_order(AccuracyOrder::First());

  /* phase change */
#if CASE < 4
  PhaseChange4 pc(mdot, mflx, q, g , f , uvw, cht,
                  time, &mixed);
  pc.set_accuracy_order(AccuracyOrder::FourthUpwind());
#else
  //PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw, time, & mixed,
  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw, time, & mixed, //  step is not used in pc
                  mixed.latent()->value(), tsat0);
#endif

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  boil::test_irun();
  boil::set_irun(1);

  int ts;
  /* load variables */
  std::vector<Scalar*> load_scalars = { &press, &c, &tpr };
  std::vector<std::string> load_scalar_names = { "press", "c", "tpr" };
  std::vector<Vector*> load_vectors = { &uvw };
  std::vector<std::string> load_vector_names = { "uvw" };

  if(boil::load_backup("time.txt",ts,time,
                       load_scalars, load_scalar_names,
                       load_vectors, load_vector_names)) {
    conc.init();
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
      c[i][j][k] = 0.0;

    const real xcent=0.0;
    const real ycent=0.0;
    const real zcent=0.0;

#if CASE < 3
    boil::setup_circle_xz(conc.color(), radius, xcent, zcent);
#elif CASE < 4
    boil::setup_sphere(conc.color(), radius, xcent, ycent, zcent);
#else
    boil::setup_sphere(c, radius, xcent, ycent, zcent);
#endif
#if CASE < 4
    conc.color().bnd_update();
    conc.color().exchange_all();
#else
    c.bnd_update();
    c.exchange_all();
#endif
#if CASE == 1
    conc.color_to_vf(conc.color(),c);
#endif
    for_vijk(c,i,j,k)
      c[i][j][k] = 1.0-c[i][j][k];

    c.bnd_update();
    c.exchange_all();
    conc.init();
    conc.totalvol();
#if CASE == 4
    set_stepfunc(c, step, sflag);
#endif

    for_vijk(c,i,j,k) {
      real dist = pow(c.xc(i)-xcent,2.0)
#if CASE > 2
                 +pow(c.yc(j)-ycent,2.0)
#endif
                 +pow(c.zc(k)-zcent,2.0);

      real radius2=70.0e-6;
      real coef1 =  2.28023719E+05;
      real coef2 = -8.71920915E+09;
      real coef3 = -7.28351952E+14;
      real coef4 =  6.46735617E+19;
      real coef5 = -1.35949950E+24;
      if (dist<=pow(radius,2.0)) {
        tpr[i][j][k] = tsat0;
      } else if(dist<=pow(radius2,2.0)) {
        real xi = sqrt(dist) - radius;
        real ttmp = coef5*pow(xi,5.0) + coef4*pow(xi,4.0) + coef3*pow(xi,3.0)
                  + coef2*pow(xi,2.0) + coef1*xi + tsat0;
        tpr[i][j][k] = std::min(tout,ttmp);
      } else {
        tpr[i][j][k] = tout;
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();
#if CASE == 4
    boil::plot->plot(uvw,c,tpr,press,mdot,step,"uvw-c-tpr-press-mdot-step",0); 
#else
    boil::plot->plot(uvw,c,tpr,press,mdot,mflx,"uvw-c-tpr-press-mdot-mflx",0); 
#endif
  }

  /* set iint */
  int iint = 1;
  boil::oout<<"iint= "<<iint<<"\n";

#ifdef FS_TEST
  #if 0
  Comp mcomp = Comp::i();
  #else
  Comp mcomp = Comp::k();
  #endif
  for_vmijk((*conc.topo.fs),mcomp,i,j,k) {
    real val = (*conc.topo.fs)[mcomp][i][j][k];
  #if 0
    real pos = (*conc.topo.fs).zc(mcomp,k);
    real sig = 1.;
  #else
    real pos  = (*conc.topo.fs).xc(mcomp,i);
    real pos2 = (*conc.topo.fs).zc(mcomp,k);
    real sig = (pos2>0.)-(pos2<0.);
  #endif
    if(boil::realistic(val)) {
      boil::oout<<i<<" "<<j<<" "<<k<<" "<<val<<" "<<sig*sqrt(radius*radius-pos*pos)
                <<" | "<<val/sig/sqrt(radius*radius-pos*pos)-1.
                <<boil::endl;
    } 
  }

  exit(0);
#endif

  MaxIter multigrid_mm_smooth1 = MaxIter(70);
  MaxIter multigrid_mm_smooth2 = MaxIter(70);
  MaxIter multigrid_mm_solve = MaxIter(110);
  MaxIter multigrid_mm_stale1 = MaxIter(-1);
  MaxIter multigrid_mm_stale2 = MaxIter(-1);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm_smooth1,multigrid_mm_smooth2,multigrid_mm_solve};
  std::array<MaxIter,3> multigrid_mstale = {multigrid_mm_stale1,multigrid_mm_stale1,multigrid_mm_stale2};

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    conc.new_time_step();

    /* store velocity & temperature */
    if(mSimple>1) {
      for_m(m)
        uvw_old(m) = uvw(m);
      tprold = tpr;
    }

    /* inner loop */
    for(int mloop=0; mloop<mSimple; mloop++) {
      /*---------------+
      |  phase change  |
      +---------------*/
      pc.update();
      ns.vol_phase_change(&f);

      /* reset body force */
      for_m(m)
        for_avmijk(xyz,m,i,j,k)
          xyz[m][i][j][k]=0.0;

      /* surface tension */
#if CASE < 4
      conc.tension(&xyz, mixed,conc.color());
#else
      conc.tension(&xyz, mixed, step);
#endif

      /*--------------------------+
      |  solve momentum equation  |
      +--------------------------*/
      /* essential for moving front */
      ns.discretize();
      pr.discretize();
      pr.coarsen();

      /* momentum */
      if(mSimple>1) {
        ns.new_time_step(uvw_old,&f);
        ns.convection();
      } else {
        ns.new_time_step(&f); //???
      }

      ns.grad(press);
      ns.solve(ResRat(1e-14));

      p = 0.0;
      //if (multigrid.wcycle(ResTol(1e-6))) OMS(converged);
      if(multigrid.cycle(Cycle::Z(),
                         Cycle::F(),
                         ResTol(-1.),
                         ResRat(1e-6),
                         multigrid_mi,
                         multigrid_mstale))
          OMS(converged);
      //pr.solve(ResRat(1e-6));
      p.exchange();
      ns.project(p);
      press += p;

      /*---------------------------+
      |  solve transport equation  |
      +---------------------------*/
#if CASE < 4
      conc.advance_with_extrapolation(conc.topo->vfold,true,ResTol(1e-9),uvw,
                                      &liquid,&uvw_1,&vapor,&uvw_2);
#else
      conc.advance();
#endif
#if CASE == 3
  #ifdef USE_MG
      conc.heaviside()->marker_gradient(*conc.topo->adens);
  #endif
#endif
#if CASE == 4
      set_stepfunc(c, step, sflag);
#endif
#if 0
    real sum = 0.0;
    int count = 0;
    for_vijk(c,i,j,k) {
      real sumplus = (*conc.topo->adens)[i][j][k]*(*conc.topo->adens).dV(i,j,k);
      if(sumplus>0.0) count++;
      sum += sumplus;
    }
    boil::cart.sum_real(&sum);
    boil::cart.sum_int(&count);
    boil::oout<<"Evaluate_adens "<<count<<" "<<sum<<" "<<4.*boil::pi*radius*radius/8.<<boil::endl;
#endif

      /*------------------------+
      |  solve energy equation  |
      +------------------------*/
#if 1
      if(mloop>0)
        tpr = tprold;
      enthFD.discretize();
      enthFD.new_time_step();
      enthFD.solve(ResRat(1e-16),"enthFD");
#else
    tprold = tpr;
    enthFD.discretize();
    //enthFD.new_time_step();
    for(int i(0); i<3;++i) {
      for(int j(0); j<3;++j) {
        enthFD.convective_time_step(tprold);
        if(j==0) {
          tprap1 = tpr;
        } else {
          tpr = (tprap1+tpr);
          tpr /= 2.;
        }
      }
      enthFD.inertial(tpr,false,Old::no);
      enthFD.solve(ResRat(1e-16));
      if(i==0) {
        tprap2 = tpr;
      } else {
        tpr = (tprap2+tpr);
        tpr /= 2.;
      }
    }
#endif
    }

    conc.totalvol();

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
    if(time.current_step() % n_per_plot == 0 || time.current_step()==1 ) {
#if CASE == 4
      boil::plot->plot(uvw,c,tpr,press,mdot,step,"uvw-c-tpr-press-mdot-step",iint); 
#else
      boil::plot->plot(uvw,c,tpr,press,mdot,mflx,"uvw-c-tpr-press-mdot-mflx",iint); 
#endif
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

    if(  boil::timer.current_min() > (wmin-12.0)
      || time.current_step()==time.total_steps()
      || time.current_time()>tend) {
      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names);
      boil::plot->plot(uvw,c,tpr,kappa,mdot,press
                      ,"uvw-c-tpr-curv-mdot-press",time.current_step());

      if(time.current_time()<tend) {
        boil::set_irun(0);
      }
      break;
    }

  } /* time loop */

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}

#if CASE == 4
/******************************************************************************/
void set_stepfunc(Scalar & c, Scalar & step, Scalar & sflag){

  const real phisurf=0.5;

  for_avijk(sflag,i,j,k) {
    sflag[i][j][k]=0;
  }
  for_vijk(c,i,j,k) {
    if(c[i][j][k]>=phisurf)
      sflag[i][j][k]=1;
  }
  /* i-direction */
  for(int i=c.si()-1; i<=c.ei(); i++){
    for_vjk(c,j,k){
       if((c[i][j][k]-phisurf)*(c[i+1][j][k]-phisurf)<=0.0){
          sflag[i  ][j][k]=2;
          sflag[i+1][j][k]=2;
       }
    }
  }
  /* j-direction */
  for(int j=c.sj()-1; j<=c.ej(); j++){
    for_vik(c,i,k){
      if((c[i][j][k]-phisurf)*(c[i][j+1][k]-phisurf)<=0.0){
          sflag[i][j  ][k]=2;
          sflag[i][j+1][k]=2;
       }
    }
  }
  /* k-direction */
  for(int k=c.sk()-1; k<=c.ek(); k++){
    for_vij(c,i,j){
       if((c[i][j][k]-phisurf)*(c[i][j][k+1]-phisurf)<=0.0){
          sflag[i][j][k  ]=2;
          sflag[i][j][k+1]=2;
       }
    }
  }
  sflag.exchange_all();
  for_avijk(c,i,j,k){
    if(sflag[i][j][k]==2){
      step[i][j][k]=c[i][j][k];
    } else {
      if(c[i][j][k]<0.5){
        step[i][j][k]=0.0;
      } else {
        step[i][j][k]=1.0;
      }
    }
  }
}
#endif
