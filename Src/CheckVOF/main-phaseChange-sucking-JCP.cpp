#include "Include/psi-boil.h"
#define USE_VOF

/* computed parameters */
const int Level = 4;  //Level = 1,2,4
const int dtLevel = 1; // 1 or 2
const int NX = 40*Level;
const int NZ = 2;

/* domain dimensions (given by problem) */
const real LX =   8.0e-3;
const real LZ =   LX*real(NZ)/real(NX);

const real tout=378.15;
const real twall=373.15;

/******************************************************************************/
main(int argc, char * argv[]) {

  // Journal of Computational Physics, 249, 127-161
  // Sucking problem

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.0*LX,1.0*LX), NX, Periodic::no() );
  Grid1D gz( Range<real>( 0.0,LZ), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gz, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::ic2());

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
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );

  /* copy b.c. from p */
  press = p.shape();
  f = p.shape();
  mdot = p.shape();
  q = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  g = c.shape();
  step = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), twall ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d);
  //vapor  .mu    (1.255e-1);
  vapor  .mu    (1.255e-5);
  vapor  .rho   (0.597);
  vapor  .cp    (2030*0.597);
  vapor  .lambda(0.025);
  //liquid.mu    (0.28e-1);
  liquid.mu    (0.28e-3);
  liquid.rho   (958.4);
  liquid.cp    (4215.9*958.4);
  liquid.lambda(0.679);
  const real latent=2258.0*1e3;
  const real tsat=373.15;

  Matter mixed(liquid, vapor, & step); //c=1: full of liquid, c=0: full of vapor

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 500*Level*dtLevel;
  //const int  ndt = 1000*Level*dtLevel;
  const int  nint = 100*Level*dtLevel;
  const real dt  = 1.0e-3/real(Level*dtLevel);

  Times time(ndt, dt);

  OPR(  NX );
  OPR(  dt );
  OPR( ndt );
  
  /*--------------------+
  |  initial condition  |
  +--------------------*/
  // read from file
  real xi[200],tprin[200];
  std::fstream input;
  input.open("input.txt");
  if( input.fail() ) {
      boil::oout<<"Error: input.txt is necesarry!!!";
  } else {
    for(int i=0; i<200; i++){
      input>>xi[i]>>tprin[i];
    }
  }
  real xint=2.0e-4;
  // color function
  for_vijk(c,i,j,k){ 
    if(c.xc(i)<xint) {
      c[i][j][k] = 0.0;
    } else {
      c[i][j][k] = 1.0;
    }
  }
  c.exchange_all();
  // temperature
  for(int i=0; i<200; i++){
    xi[i]+=xint;
  }
  real xiend=xi[199];
  for_vijk(tpr,i,j,k){
    real xtmp=tpr.xc(i);
    if(xtmp>=xiend){
      tpr[i][j][k]=tout;
    } else if (xtmp<=xint) {
      tpr[i][j][k]=tsat;
    } else {
      for(int ii=0; ii<199; ii++){
        if((xtmp-xi[ii])*(xtmp-xi[ii+1])<=0) {
          real we1=(xi[ii+1]-xtmp)/(xi[ii+1]-xi[ii]);
          real we2=(xtmp-xi[ii])/(xi[ii+1]-xi[ii]);
          tpr[i][j][k]=we1*tprin[ii]+we2*tprin[ii+1];
        }
      }
    }
  }
  tpr.exchange_all();

  boil::plot->plot(uvw, c, tpr, mdot, "uvw-c-tpr-mdot", 0);

  for_avijk(c,i,j,k){
    step[i][j][k]=c[i][j][k];
  }

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set(TimeScheme::backward_euler());
  Pressure pr(p, f, uvw, time, solver, &mixed);
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, &mixed ,tsat);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
#ifdef USE_VOF
  VOF conc  (c,  g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
#endif
  conc.front_minmax();
  conc.totalvol();
  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw, time, &mixed , latent, tsat);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  /*--------+
  |  model  |
  +--------*/
  Model tm;

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

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    cold=c;
    conc.advance();
#if 1
    for_vi(c,i) {
      real csum=0.0;
      real isum=0;
      for_vjk(c,j,k) {
        csum += c[i][j][k];
        isum++;
      }
      for_vjk(c,j,k) {
        c[i][j][k]=csum/real(isum);
      }
    }
    c.bnd_update();
    c.exchange_all();
#endif
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
    enthFD.solve(ResRat(1e-16),"enthalpy");

#if 1
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
    ns.solve(ResRat(0.01));

    p = 0.0;
    //if (multigrid.vcycle(1e-8)) OMS(converged);
    if (multigrid.vcycle(ResRat(1e-4))) OMS(converged);
    //pr.solve(1e-12);
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
#endif
    if(time.current_step() % (nint)==0 ||time.current_step()==1) {
      boil::plot->plot(uvw,c,tpr,mdot,"uvw-c-tpr-mdot",time.current_step());
    }
  }

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

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-JCP-phaseChange-sucking.cpp,v 1.3 2018/04/30 08:45:18 sato Exp $'/
+-----------------------------------------------------------------------------*/
