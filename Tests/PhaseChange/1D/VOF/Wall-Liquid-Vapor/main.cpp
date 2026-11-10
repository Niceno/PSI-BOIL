#include "Include/psi-boil.h"

/******************************************************************************/
int main(int argc, char * argv[]) {

  boil::timer.start();

  /* computed parameters */
  const int Level = 2;    // 1
  const int NX = 16*Level;

  /* domain dimensions (given by problem) */
  const real LX = 0.001; //Hardt
  const real DX = LX/real(NX);
  const real thick = LX*0.5;  // film thicknes

  const real Twall=200;
  const real Tsat=100;

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.0*LX,1.0*LX), NX, Periodic::no() );
  Grid1D gz(DX);

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gz, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

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
  Vector uvw(d), xyz(d), uvw_1(d), uvw_2(d); // vel
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d), cold(d); // concentration
  Scalar tpr(d), q  (d); // temperature
  Scalar mdot(d), mflx(d);        // phase-change

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
    uvw_1(m) = uvw(m).shape();
    uvw_2(m) = uvw(m).shape();
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );

  /* copy b.c. from p */
  press = p.shape();
  f = p.shape();
  mdot = p.shape();
  mflx = p.shape();
  q = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );

  g = c.shape();
  cold = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), Twall ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), Tsat ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d);
  vapor  .mu    (1.255e-5);
  vapor  .rho   (0.597);
  vapor  .cp    (2030*0.597);
  vapor  .lambda(0.025);
  liquid.mu    (0.28e-3);
  liquid.rho   (958.4);
  liquid.cp    (4215.9*958.4);
  liquid.lambda(0.679);

  Matter mixed(liquid, vapor,& c); //c=1: full of liquid, c=0: full of vapor
  mixed.sigma(2.3610e-2);
  mixed.latent(2258.0*1e3);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 100000*Level;
  const real dt  = 1.0e-4/real(Level);
  const int  nint = 10;
  const real tint = 2.0e-2;    // time interval to output
  const real cfl_limit = 0.25; // maximum Courant number
  Times time(ndt, dt);
  boil::oout<<"dt= "<<dt<<"\n";

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  //boil::test_irun();
  //boil::set_irun(1);
 
  /*--------------------+
  |  initial condition  |
  +--------------------*/
  // VOF
  for_vijk(c,i,j,k){
    if(c.xn(i+1)<=thick) {
      c[i][j][k] = 1.0;
    } else if(c.xn(i)<thick && thick<c.xn(i+1)) {
      c[i][j][k] = (thick-c.xn(i))/c.dxc(i);
    } else {
      c[i][j][k] = 0.0;
    }
  }
  c.bnd_update();
  c.exchange_all();

  // temperature
  for_vijk(tpr,i,j,k) {
    if(c.xc(i) < thick) {
      tpr[i][j][k] = Twall + (Tsat-Twall)/thick * c.xc(i);
    } else {
      tpr[i][j][k] = Tsat;
    }
  }
  tpr.bnd_update();
  tpr.exchange_all();

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.convection_set(TimeScheme::forward_euler());
  ns.diffusion_set(TimeScheme::backward_euler());
  Pressure pr(p, f, uvw, time, solver, &mixed);
  VOF conc  (c,  g, kappa, uvw_1, time, solver);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(1);

  conc.init();
  c.exchange_all();
  conc.ancillary();
  conc.front_minmax();
  conc.totalvol();

  TIF tsat(Tsat);
  CommonHeatTransfer cht(tpr,conc.topo,tsat,&mixed);
  EnthalpyFD enthFD(tpr, q, uvw, uvw_1, uvw_2, time, solver, &mixed, cht);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  //enthFD.set_no_solid_acceleration(false);
  enthFD.set_flux_accuracy_order(AccuracyOrder::First());
  PhaseChange4 pc(mdot, mflx, q, g , f , uvw, cht,
                  time, &mixed);
  pc.set_accuracy_order(AccuracyOrder::FourthUpwind());
  pc.set_unconditional_extrapolation(false);
  pc.set_discard_points_near_interface(false);

  // output initial condition
  boil::plot->plot(uvw, c, tpr, mdot, "uvw-c-tpr-mdot", 0);

  /* set iint (counter for TEC file) */
  int iint = int(time.current_time()/tint) + 1;
  boil::oout<<"iint= "<<iint<<"\n";

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {
    /* reset RHS of momentum equation */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
	xyz[m][i][j][k]=0.0;

    /* new time step */
    conc.new_time_step();

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-14));

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();                // calculate mdot
    ns.vol_phase_change(&f);    // set RHS due to mdot

    /* essential for moving front */
    ns.discretize();
    pr.discretize();

    /* momentum */
    ns.new_time_step();
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

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    cold=c;
    conc.advance_with_extrapolation(true,ResTol(1e-7),uvw,f,
                                    &liquid,&uvw_1,&vapor,&uvw_2);
    conc.front_minmax();
    conc.totalvol();

    /*-------------+
    |  dt control  |
    +-------------*/
    real max_cfl = ns.cfl_max();
    boil::oout<<"max_cfl= "<<time.current_time()<<" "<<max_cfl<<"\n";
    time.control_dt(max_cfl, cfl_limit, dt);

    /*---------------+
    |  post-process  |
    +---------------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,c,tpr,mdot,"uvw-c-tpr-mdot",iint);
      iint++;
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::oout <<"$ ./Boil > log.txt\n";
  boil::oout <<"$ cat log.txt |grep xmin > front.out\n";
  boil::oout <<"$ gnuplot\n";
  boil::oout <<"  plot sqrt(2.5e-7-6.6917e-8*x),\"front.out\" u 2:3 w l\n";

  boil::timer.stop();
  boil::timer.report();

}

