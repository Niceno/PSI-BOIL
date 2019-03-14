#include "Include/psi-boil.h"
#define USE_VOF

const int Level=2;

/* computed parameters */
const int NX = 50*Level;
const int NZ = 2;

/* domain dimensions (given by problem) */
const real LX =   0.001; //Hardt
const real LZ =   LX*real(NZ)/real(NX);

const real Tout=110;
const real Tsat=100;

/******************************************************************************/
main(int argc, char * argv[]) {

  // Journal of Computational Physics, 249, 127-161
  // Stefan problem

  boil::timer.start();

  boil::oout<<"Edit cipcsl2_sharpen.cpp!\n";
  boil::oout<<"#if 1 (for 1D)\n";

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.0*LX,1.0*LX), NX, Periodic::no() );
  Grid1D gz( Range<real>( 0.0,LZ), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gz, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // pressure
  Scalar c  (d), g  (d), kappa(d), cold(d); // concentration
  Scalar press(d);
  Scalar tpr(d), q  (d), step(d); // temperature
  Scalar mdot(d);        // phase-change
  Scalar nx(d),ny(d),nz(d);        // phase-change

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  /* copy b.c. from p */
  press = p.shape();
  f = p.shape();
  mdot = p.shape();
  q = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  g = c.shape();
  step = c.shape();
  cold = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), Tout ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), Tsat ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  OPR( Tout );

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
  const real latent=2258.0*1e3;

  Matter mixed(liquid, vapor,& step); //c=1: full of liquid, c=0: full of vapor
  mixed.sigma(2.3610e-2);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 5000*Level;
  const int  nint = 500;
  const real dt  = 2.0e-5/real(Level);
  Times time(ndt, dt);
  time.print_time(false);

  OPR(  NX );
  OPR(  dt );
  OPR( ndt );
  
  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

  for_vijk(c,i,j,k) {
    if(c.xc(i)<LX/real(NX)) {
      c[i][j][k]=0.0;
    } else if (i==2) {
      c[i][j][k]=0.99;
    } else {
      c[i][j][k]=1.0;
    }
  }
  c.exchange_all();

  for_avijk(tpr,i,j,k){
    tpr[i][j][k] = c[i][j][k]*Tsat
                 + (1.0-c[i][j][k])*Tout;
  }
  boil::plot->plot(uvw, c, tpr, mdot, "uvw-c-tpr-mdot", 0);

  for_avijk(c,i,j,k){
    step[i][j][k]=c[i][j][k];
  }

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  Pressure pr(p, f, uvw, time, solver, &mixed);
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, &mixed, Tsat);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
#ifdef USE_VOF
  VOF conc (c,  g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
#endif
  conc.front_minmax();
  conc.totalvol();
  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw, time, &mixed, latent, Tsat);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);
  multigrid.min_cycles(3);

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "########################" << boil::endl;

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();
    ns.vol_phase_change(&f);

    /* essential for moving front */
    ns.discretize();
    pr.discretize();

    /* momentum */
    ns.cfl_max();
    ns.new_time_step();

    for_m(m)
      for_avmijk(xyz,m,i,j,k)
	xyz[m][i][j][k]=0.0;

    ns.grad(press);
    ns.solve(ResRat(0.001));
    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-8))) OMS(converged);
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();

#if 1
    /* this part does not work in parallel computation */
    real dltp = press[1][1][1];
    for_vijk(press,i,j,k)
      press[i][j][k] -= dltp;
#endif

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    cold=c;
    conc.advance();
    conc.front_minmax();
    boil::oout<<"x-min-front= "<<time.current_time()<<" "
              <<conc.get_xminft()<<"\n";
    conc.totalvol();

    for_avijk(c,i,j,k){
      step[i][j][k]=c[i][j][k];
    }

    pc.modify_vel(uvw,c,cold);

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-4));

    std::cout<<"c[2.5]= "<<c[2][1][1]<<" "<<c[3][1][1]<<" "
             <<0.5*(c[2][1][1]+c[3][1][1])<<" "
             <<uvw[Comp::u()][3][1][1]<<"\n";

    if(time.current_step() % (nint)==0 ||time.current_step()==1) {
      boil::plot->plot(uvw,c,tpr,mdot,"uvw-c-tpr-mdot",time.current_step());
    }
    if(time.current_step() == 500) {

      std::ofstream fout;
      fout.open("velocity-u.txt");
      Comp m=Comp::u();
      for_vmi(uvw,m,i) {
        fout << uvw.xc(m,i) << "  " << uvw[m][i][1][1] << "\n";
      }
      fout.close();

      fout.open("color.txt");
      for_vi(c,i) {
        fout << c.xc(i) << "  " << c[i][1][1] << "\n";
      }
      fout.close();

      fout.open("mdot.txt");
      for_vi(c,i) {
        fout << c.xc(i) << "  " << mdot[i][1][1] << "\n";
      }
      fout.close();

      fout.open("temperature.txt");
      for_vi(tpr,i) {
        fout << tpr.xc(i) << "  " << tpr[i][1][1] << "\n";
      }
      fout.close();
    }
  }

#if 1
  std::ofstream fout;
  fout.open("velocity-u.txt");
  Comp m=Comp::u();
  for_vmi(uvw,m,i) {
       fout << uvw.xc(m,i) << "  " << uvw[m][i][1][1] << "\n";
  }
  fout.close();

  fout.open("color.txt");
  for_vi(c,i) {
       fout << c.xc(i) << "  " << c[i][1][1] << "\n";
  }
  fout.close();

  fout.open("temperature.txt");
  for_vi(tpr,i) {
       fout << tpr.xc(i) << "  " << tpr[i][1][1] << "\n";
  }
  fout.close();

  fout.open("mdot.txt");
  for_vi(tpr,i) {
       fout << mdot.xc(i) << "  " << mdot[i][1][1] << "\n";
  }
  fout.close();
#endif

  boil::oout << "finished" << boil::endl;

  boil::oout <<"#  cat log.txt |grep x-min-front > front.out\n";
  boil::oout <<"#  Use gnuplot\n";
  boil::oout <<"#  plot 2*6.695e-2*sqrt(0.025*x/(0.597*2030)) ,\"front.out\" u 2:3 w l\n";

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-JCP-phaseChange-stefan.cpp,v 1.3 2018/04/30 08:45:18 sato Exp $'/
+-----------------------------------------------------------------------------*/
