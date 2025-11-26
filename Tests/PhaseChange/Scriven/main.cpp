#include "Include/psi-boil.h"
#define USE_VOF

/* computed parameters */
const int gLevel = 2;  //grid level=2,4,8
const int NX = 24*gLevel;

/* domain dimensions (given by problem) */
const real LX = 187.5e-6;
const real radius=50.0e-6;

const real tsat0=0.0;
const real tout=tsat0+1.25;
const real dtsatdp=3.0/4.0;
const real pi=acos(-1.0);
real tsat = tsat0;

void set_stepfunc(Scalar & sca, Scalar & scb, Scalar & scc);
real interpT(const std::vector<real>& r,
             const std::vector<real>& T, real rp);


/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc==1){
    boil::oout<<"An argument is required!"<<"\n";
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
  Grid1D gx( Range<real>(0.0,LX), NX, Periodic::no(),
             BndGrid::symmetry(), BndGrid::extrapolate());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gx);

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
  Scalar tpr(d), q  (d); // temperature
  Scalar step(d), sflag(d); // heaviside function
  Scalar mdot(d);        // phase-change

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /* copy b.c. from p */
  press = p.shape();
  f     = p.shape();
  mdot  = p.shape();
  q     = p.shape();
  kappa = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  g = c.shape();
  step = c.shape();
  sflag = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), tout ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), tout ) );

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
#ifdef USE_VOF
  Matter mixed(liquid, vapor, & c);
#else
  Matter mixed(liquid, vapor, & step);
#endif
  mixed.sigma(5.9e-2);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt =  5000*gLevel;
  const int  nint =  350*gLevel;
  const real tint = 2.0e-5;
  const real dxmin = LX/real(NX);
  const real dt = 100.0 
                * pow(0.5*pow(dxmin,3.0)/(2.0*3.1415*mixed.sigma()->value()), 0.5);
  Times time(ndt, dt);
  time.end(0.0018);
  boil::oout<<"main:dt= "<<dt<<"\n";
  time.print_time(false);
  
  /*-----------------+
  |  define solvers  |
  +-----------------*/
#ifdef USE_VOF
  VOF conc  (c,  g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_itsharpen(10);
#endif

  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw, time, &mixed, latent, tsat);

  Momentum ns( uvw, xyz, time, solver, &mixed);

  Pressure pr(p, f, uvw, time, solver, &mixed);

  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, &mixed, tsat);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

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

  int ts;
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
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
  } else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    /*------------------------------+
    |  initial temperature profile  |
    +------------------------------*/
    const std::string filename = "dT1.25_R50microns.txt";

    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "Error: cannot open file " << filename << "\n";
        return 1;
    }

    std::vector<real> r, T;
    std::string line;

    while (std::getline(fin, line)) {
      // skip empty lines
      if (line.empty()) continue;

      // skip header or any non-numeric-leading line
      // (allows for "r   T" at top)
      char c0 = line.find_first_not_of(" \t") != std::string::npos
                  ? line[line.find_first_not_of(" \t")]
                  : '\0';
      if (!(std::isdigit(c0) || c0=='-' || c0=='+' || c0=='.')) continue;

      std::istringstream iss(line);
      real rv, Tv;
      if (iss >> rv >> Tv) {
          r.push_back(rv);
          T.push_back(Tv-100.0);
      }
      // else: silently ignore malformed lines
    }

    fin.close();

    std::cout << "Read " << r.size() << " rows.\n";
    if (!r.empty()) {
        boil::oout << "First row: r=" << r.front() << ", T=" << T.front() << "\n";
        boil::oout << "Last row : r=" << r.back()  << ", T=" << T.back()  << "\n";
    }

    /*--------------------+
    |  initial condition  |
    +--------------------*/
    for_vijk(c,i,j,k) 
      c[i][j][k] = 0.0;

    const real xcent=0.0;
    const real ycent=0.0;
    const real zcent=0.0;
    for_vijk(c,i,j,k) {
      real dist = sqrt( pow(c.xc(i)-xcent,2.0)
                       +pow(c.yc(j)-ycent,2.0)
                       +pow(c.zc(k)-zcent,2.0));

      if (dist<=radius) {
        tpr[i][j][k] = tsat;
      } else if(dist<=r[r.size()-1]) {
        tpr[i][j][k] = interpT(r, T, dist);
      } else {
        tpr[i][j][k] = tout;
      }
    }

#ifdef USE_VOF
    boil::setup_sphere(c, radius, xcent, ycent, zcent, 20);
#else
    for_vijk(c,i,j,k) {
      real dist = sqrt( pow(c.xc(i)-xcent,2.0)
                       +pow(c.yc(j)-ycent,2.0)
                       +pow(c.zc(k)-zcent,2.0));
      real dd = -dist+radius;
      real eps = dxmin*1.5;
      if(dd<-eps){
        c[i][j][k]=0.0;
      } else if (dd<eps){
        c[i][j][k]=0.5+dd/(2.0*eps)+1/(2*pi)*sin(pi*dd/eps);
      } else {
        c[i][j][k]=1.0;
      }
    }
#endif
    for_vijk(c,i,j,k)
      c[i][j][k] = 1.0-c[i][j][k];

    c.exchange_all();
    tpr.exchange_all();

    conc.init();

#ifndef USE_VOF
#if 1
    for_avijk(c,i,j,k){
      step[i][j][k]=c[i][j][k];
    }
#else
  set_stepfunc(c, step, sflag);
#endif
#endif

    boil::plot->plot(uvw,c,tpr,press,mdot,"uvw-c-tpr-press-mdot",0);

  }
  input.close();

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;
  boil::oout<<"iint= "<<iint<<"\n";

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "########################" << boil::endl;
    boil::oout << " dt= " << time.dt() << boil::endl;

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();
    ns.vol_phase_change(&f);

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthFD");

    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    cold = c;
    conc.advance();
    conc.tension(&xyz, mixed);
    conc.front_minmax();
    boil::oout<<"x-min-front= "<<time.current_time()<<" "<<conc.get_xminft()
              <<" "<<conc.get_xmaxft()<<"\n";
    conc.totalvol();

#if 1
  for_avijk(c,i,j,k){
    step[i][j][k]=c[i][j][k];
  }
#else
  set_stepfunc(c, step, sflag);
#endif

    pc.modify_vel(uvw,c,cold);

    /* essential for moving front */
    ns.discretize();
    pr.discretize();

    /* momentum */
    ns.cfl_max();
    ns.new_time_step();

    ns.grad(press);
    ns.solve(ResRat(1.0e-14));
    p = 0.0;
    boil::oout<<"main:multigrid\n";
    multigrid.vcycle(ResRat(1e-6));
    boil::oout<<"main:p.exchange\n";
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(pmin>press[i][j][k]) pmin=press[i][j][k];
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k)
      press[i][j][k] -= pmin;

    /* dt control */
    time.control_dt(ns.cfl_max(),0.10,dt);
    boil::oout<<"main:after dt control\n";

    /* output data */
    if((time.current_time()) / (tint) >= real(iint) ) {
      tpr.exchange_all();
      boil::plot->plot(uvw,c,tpr,press,mdot,"uvw-c-tpr-press-mdot",iint);
      iint++;
    }

    if((time.current_step()) % (nint)==0 ) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
    }

    if( boil::timer.current_min() > wmin ) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
      std::fstream output;
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      output.open("run.txt", std::ios::out);
      output << 0 << boil::endl;
      output.close();
      boil::timer.stop();
      boil::timer.report();
      uvw  .rm("uvw", ts);
      press.rm("press", ts);
      //conc .rm("conc", ts);
      tpr  .rm("tpr", ts);
      exit(0); 
    }
    boil::oout<<"main:End of Time Loop\n";
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}

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

real interpT(const std::vector<real>& r,
               const std::vector<real>& T,
               real rp)
{
    if (r.size() != T.size() || r.empty())
        throw std::runtime_error("interpT: r and T must be same nonzero size.");

    // range check (as you requested)
    const real rmin = r[0], rmax = r[r.size()-1];
    if (rp < rmin ) return T[0];
    if (rp > rmax ) return T[T.size()-1];
    if (rp < rmin || rp > rmax)
        throw std::out_of_range("interpT: rp is outside [5e-5, 7e-5].");

    // find first element >= rp
    auto it = std::lower_bound(r.begin(), r.end(), rp);

    // exact match or boundary
    if (it == r.begin()) return T.front();
    if (it == r.end())   return T.back();
    size_t i = static_cast<size_t>(it - r.begin());
    if (r[i] == rp) return T[i];

    // linear interpolation between i-1 and i
    real r1 = r[i - 1], r2 = r[i];
    real T1 = T[i - 1], T2 = T[i];

    return T1 + (T2 - T1) * (rp - r1) / (r2 - r1);
}
