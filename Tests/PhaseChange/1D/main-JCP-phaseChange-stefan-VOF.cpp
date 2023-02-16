#include "Include/psi-boil.h"

#define USE_VOF

/******************************************************************************/
int main(int argc, char * argv[]) {

  boil::timer.start();

  if(argc<4) {
    boil::oout<<"Three command line arguments required!"<<"\n";
    boil::oout<<"./Boil wmin level cfl (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  const int Level=atoi(argv[2]);
  boil::oout<<"level= "<<Level<<"\n";

  const real cfl_limit=atof(argv[3]);
  boil::oout<<"cfllim= "<<cfl_limit<<"\n";

#ifdef USE_VOF
  boil::oout<<"#Stephan problem using VOF\n";
#else
  boil::oout<<"#Stephan problem using CIPCSL2\n";
#endif

  /* computed parameters */
  const int NX = 50*Level;
  const int mSimple = 1;

  /* dt settings */
  const real initdtcoef = 1./10.;

  /* domain dimensions (given by problem) */
  const real LX = 0.001; //Hardt
  const real DX = LX/real(NX);

  const real betasol = 0.06694560396502779;
  const real t0 = 30e-3;
  const real xint = 105.32888985625091e-6;
  const real x0 = 0.0;//xint-5.01*DX*Level;

  const real Twall=110;
  const real Tsat=100;

#ifndef USE_VOF
  boil::oout<<"Edit cipcsl2_sharpen.cpp!\n";
  boil::oout<<"#if 1 (for 1D)\n";
  const int NZ = 3;
  const real LZ = DX*real(NZ);
#endif

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();
  //boil::plot = new PlotTEC(AsNodes::no(), Buffers::yes());

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(0.0*LX+x0,1.0*LX+x0), NX, Periodic::no() );
#ifndef USE_VOF
  Grid1D gz( Range<real>( 0.0,LZ), NZ, Periodic::yes() );
#else
  Grid1D gz(DX);
#endif

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gz, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  //Krylov * solver = new CG(d, Prec::ic2());
  Krylov * solver = new BiCGS(d, Prec::di());
  //Linear * solver = new Jacobi(d);

  MaxIter multigrid_mm_smooth1 = MaxIter(110);
  MaxIter multigrid_mm_smooth2 = MaxIter(110);
  MaxIter multigrid_mm_solve = MaxIter(110);
  MaxIter multigrid_mm_stale1 = MaxIter(-1);
  MaxIter multigrid_mm_stale2 = MaxIter(-1);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm_smooth1,multigrid_mm_smooth2,multigrid_mm_solve};
  std::array<MaxIter,3> multigrid_mstale = {multigrid_mm_stale1,multigrid_mm_stale1,multigrid_mm_stale2};

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d), uvw_1(d), uvw_2(d), uvw_old(d); // vel
  Scalar p  (d), f  (d); // pressure
  Scalar c  (d), g  (d), kappa(d), cold(d); // concentration
  Scalar press(d);
  Scalar tpr(d), q  (d), tprold(d); // temperature
  Scalar mdot(d);        // phase-change
  Scalar mflx(d);        // phase-change
  Scalar sflag(d);        // phase-change

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
#ifndef USE_VOF
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
#else
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
#endif
    uvw_1(m) = uvw(m).shape();
    uvw_2(m) = uvw(m).shape();
    uvw_old(m) = uvw(m).shape();
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
#ifndef USE_VOF
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
#else
  p.bc().add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
#endif

  /* copy b.c. from p */
  press = p.shape();
  f = p.shape();
  mdot = p.shape();
  mflx = p.shape();
  q = p.shape();

  c = p.shape();
  g = c.shape();
  cold = c.shape();
  sflag = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), Twall ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), Tsat ) );
#ifndef USE_VOF
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
#else
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
#endif
  tprold = tpr.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  const real muv = 1.255e-5;
  const real rhov = 0.597;
  const real lambdav = 0.025;
  const real cpv = 2030.*rhov;
  const real alpv = lambdav/cpv;

  const real mul = 0.28e-3;
  const real rhol = 958.4;
  const real lambdal = 0.679;
  const real cpl = 4215.9*rhol;
  const real alpl = lambdal/cpl;

  Matter vapor(d), liquid(d);
  vapor  .mu    (muv);
  vapor  .rho   (rhov);
  vapor  .cp    (cpv);
  vapor  .lambda(lambdav);
  liquid.mu    (mul);
  liquid.rho   (rhol);
  liquid.cp    (cpl);
  liquid.lambda(lambdal);

  Matter mixed(liquid, vapor,& c); //c=1: full of liquid, c=0: full of vapor
  mixed.sigma(2.3610e-2);
  mixed.latent(2258.0*1e3);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 50000*Level;
  const int  nint = 100;
  const real dt  = cfl_limit/9./real(Level);

  Times time(ndt, dt);
  time.set_coef_dec(0.75);
  time.set_dt(dt*initdtcoef);

  OPR(  NX );
  OPR(  LX );
  OPR(  dt );
  OPR( ndt );

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  boil::test_irun();
  boil::set_irun(1);
 
  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vijk(c,i,j,k)
    c[i][j][k] = 1.0;
  for_vijk(c,i,j,k) {
    if(c.xn(i+1)<xint) {
      c[i][j][k]=0.0;
    } else if(c.xn(i)<xint) {
      c[i][j][k]=1.-(xint-c.xn(i))/c.dxc(i);
    } else {
      c[i][j][k] = 1.0;
    }
  }
  c.bnd_update();
  c.exchange_all();

  /* temperature */
  const real gamma = rhov*sqrt(alpv)/rhol/sqrt(alpl);
  auto tsol = [&](const real x, const real t) {
    return Twall + (Tsat-Twall)/erf(betasol)*erf(x/2./sqrt(alpv*t));
  };
  for_vijk(tpr,i,j,k) {
    real xtmp=tpr.xc(i);
    if(xtmp>=xint) {
      tpr[i][j][k] = Tsat;
    } else {
      tpr[i][j][k] = tsol(xtmp,t0);
    }
  }
  tpr.bnd_update();
  tpr.exchange_all();

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  if(mSimple>1)
    ns.convection_set(TimeScheme::backward_euler()); //ns.convection is mandatory
  else
    ns.convection_set(TimeScheme::forward_euler());
  ns.diffusion_set(TimeScheme::backward_euler());
  Pressure pr(p, f, uvw, time, solver, &mixed);
#ifndef USE_VOF
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
#else 
  VOF conc  (c,  g, kappa, uvw_1, time, solver);
#endif
  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(1);

  conc.init();
  c.exchange_all();
  conc.ancillary();
  //conc.front_minmax();
  //conc.totalvol();

  TIF tsat(Tsat);
  CommonHeatTransfer cht(tpr,conc.topo,tsat,&mixed);
  EnthalpyFD enthFD(tpr, q, uvw, uvw_1, uvw_2, time, solver, &mixed, cht);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  enthFD.set_flux_accuracy_order(AccuracyOrder::First());
  PhaseChange4 pc(mdot, mflx, q, g , f , uvw, cht,
                  time, &mixed);
  pc.set_accuracy_order(AccuracyOrder::FourthUpwind());
  //pc.set_accuracy_order(AccuracyOrder::Second());
  pc.set_unconditional_extrapolation(false);
  pc.set_discard_points_near_interface(false);

  boil::plot->plot(uvw, c, tpr, mdot, "uvw-c-tpr-mdot", 0);

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /* store velocity & temperature */
    if(mSimple>1) {
      for_m(m)
        uvw_old(m) = uvw(m);
      tprold = tpr;
    }

    /* new time step */
    conc.new_time_step();

    /* inner loop */
    for(int mloop=0; mloop<mSimple; mloop++) {
      /*---------------+
      |  phase change  |
      +---------------*/
      pc.update();

      ns.vol_phase_change(&f);
      real massflux_heat = pc.get_smdot();
      massflux_heat /= conc.topo->get_totarea();
      boil::oout<<"mflux= "<<time.current_time()<<" "
                           <<massflux_heat<<" "
                           <<boil::endl;


#if 0
      /* essential for moving front */
      ns.discretize();
      pr.discretize();

      /* momentum */
      if(mSimple>1) {
        ns.new_time_step(uvw_old,&f);
        ns.convection();
      } else {
        ns.new_time_step(&f);
      }

      for_m(m)
        for_avmijk(xyz,m,i,j,k)
          xyz[m][i][j][k]=0.0;

      ns.grad(press);
      ns.solve(ResRat(1e-14));

      p = 0.0;
      //if (multigrid.fcycle(ResTol(1e-7))) OMS(converged);
      if(multigrid.cycle(Cycle::Z(),
                         Cycle::F(),
                         ResTol(1e-7),
                         ResRat(-1.),
                         multigrid_mi,
                         multigrid_mstale))
          OMS(converged);

      p.exchange();
      ns.project(p);
      press += p;
      press.exchange();
#else
      Comp m = Comp::u();
      real xpos = conc.topo->get_xmaxft();
      for_vmijk(uvw,m,i,j,k) {
        if(uvw.xc(m,i)>xpos) {
          uvw[m][i][j][k] = massflux_heat*(1./rhov-1./rhol);
        } else {
          uvw[m][i][j][k] = 0.0;
        }
      }
      uvw.bnd_update_nooutlet();
      ns.outlet();
      uvw.exchange_all();
      for_avmijk(uvw,m,i,j,k) {
        uvw_1[m][i][j][k] = massflux_heat*(1./rhov-1./rhol);
        uvw_2[m][i][j][k] = 0.0;
      }
#endif

      /*---------------------------+
      |  solve transport equation  |
      +---------------------------*/
#ifndef USE_VOF
      conc.advance();
#else
      //conc.advance_with_extrapolation(conc.topo->vfold,true,ResTol(1e-7),uvw,
      //                                &liquid,&uvw_1,&vapor,&uvw_2);
      conc.advance(conc.topo->vfold);
#endif

#ifndef USE_VOF
      //step = c;
      //pc.modify_vel(uvw,c,conc.topo->vfold);
#endif

      /*------------------------+
      |  solve energy equation  |
      +------------------------*/
      if(mloop>0)
        tpr = tprold;
      enthFD.discretize();
      enthFD.new_time_step();
      enthFD.solve(ResRat(1e-10),"enthFD");
    }

    conc.front_minmax();
    conc.totalvol();
    time.control_dt(ns.cfl_max(),cfl_limit,dt);

    if(time.current_step() % ndt == 0 ||time.current_step()==1) {
      boil::plot->plot(uvw,c,tpr,mdot,"uvw-c-tpr-mdot",time.current_step());
    }

    if(conc.get_xmaxft()>500e-6) {
      boil::plot->plot(uvw,c,tpr,mdot,"uvw-c-tpr-mdot",time.current_step());
      break;
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::oout <<"#  Use gnuplot\n";
  boil::oout <<"#  plot 2*6.695e-2*sqrt(0.025*x/(0.597*2030)) ,\"front.out\" u 2:3 w l\n";

  boil::timer.stop();
  boil::timer.report();

}
