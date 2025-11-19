#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#define OPT_RACK
using namespace std;

#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

const int  gLevel = 16; //32

/* Domain */
const int  NX1 = 16*gLevel;
const int  NX2 =  2*gLevel;
const int  NY  =  8*gLevel;
const int  NZ  =  8*gLevel;
const real LY  =  0.005;
const real LX1 =  0.01;
const real LX3 =  0.012;
const real LZ  =  0.00589;

/* parameter for boundary conditions */
const real Uc = 1.127234319;  // average velocity (m/s)

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
  const real dx = LX1/real(NX1);
  Grid1D gx1( Range<real>(0.0, LX1), NX1, Periodic::no());
  Grid1D gx2( Range<real>(LX1, LX3), Range<real>(1.2*dx,2*dx), NX2, Periodic::no());
  Grid1D gx( gx1, gx2, Periodic::no());
  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), NY, Periodic::yes());

  Grid1D gz1( Range<real>(0.0,0.5*LZ), Range<real>(dx/4.0,dx),NZ*3/4, Periodic::no());
  Grid1D gz2( Range<real>(0.5*LZ,LZ),  Range<real>(1.2*dx,4.0*dx), NZ/4, Periodic::no());
  Grid1D gz( gz1, gz2, Periodic::no(), BndGrid::wall(), BndGrid::symmetry());

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Domain dom(gx, gy, gz, & floor);
  //boil::plot->plot(dom,"dom");
  dom.save("grid16.grd");

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(dom), xyz(dom);             // velocity
  Scalar press(dom), p  (dom), f  (dom); // pressure
  Scalar wd(dom), mu_t(dom), yplus(dom);// eddy viscosity
  Scalar sdummy(dom);  // temporary or nousage
  real *** copy_planeVec;
  //alloc3d( &copy_planeVec, 3, uvw.nj()+1, uvw.nk()+1);
  //real copy_planeVec[3][uvw.nj()+1][uvw.nk()+1];

#ifdef OPT_RACK
  #include "monitor_define.cpp"
#endif

  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::insert() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  }
  press.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  f    = p.shape();
  yplus = p.shape();
  mu_t  = p.shape();
  wd    = p.shape();
  sdummy  = p.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter liquid(dom); // zircaloy(dom)
  liquid.mu       (1.0e-6);
  liquid.rho      (1.0e+3);
  real const mu_t_max =100.0*liquid.mu()->value();  // previously 10.0

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt  = 100000;
  const real tint = 1.0e-3;  // time interval for dat-output (s)
  const real tint2 = 2.0e-5; // time interval for monitoring rack (s)
  const int  nint = 5000;    // interval for bck-output
  const real dxmin = dom.dxyz_min(); // minimum size of cell in liquid
  const real dt_max = 1e-4;   // max dt (s)
  const real dt0 = 1e-5;      // initial dt (s)
  const real cfl_limit=0.25;  // CFL limit
  Times time(ndt, dt0);
  time.print_time(false);
  time.set_coef_dec(0.5);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(dom, Prec::ic2());

  Pressure pr( p,   f,   uvw, time, solver, &liquid );
  Momentum ns( uvw, xyz,      time, solver, &liquid );
  ns.convection_set(TimeScheme::adams_bashforth());
  ns.convection_set(ConvScheme::central());
  ns.diffusion_set(TimeScheme::crank_nicolson());

  AC sol( &pr );
  sol.stop_if_diverging(true);
  sol.min_cycles(3);
  sol.max_cycles(10);

  /* wall distance */
  Distance di(wd, sdummy, uvw, time, solver);

  /* eddy viscosity */
  Model tm;

  /* SEM Inlet */
  SEMInlet sem;
  sem.set_func("tent");  // another possibility:"gauss", but need to debug it.

  /*-----------------+
  |  inlet velocity  |
  +-----------------*/
  // read fluent-2DR3.txt
  #include "read_fluent.cpp"

  //SEM: sig_rep: representative sig
  real sum_sig=0.0;
  int n_sig=0;
  real sig_min=1e+10;
  for (int KC=0; KC<dom.gik(); KC++){
    real ztmp = dom.zc_global(KC+boil::BW);
    if(ztmp>0.0){
      sum_sig+=sig_global[KC];
      n_sig++;
      if (sig_min>sig_global[KC]) sig_min=sig_global[KC];
    }
  }
  const real sig_rep = sum_sig/real(n_sig);  // average sigma

  //SEM setup
  const real alpha_SEM = 1.0;
  sem.setup(-0.5*LY,0.5*LY,0.0,LZ,sig_rep,Uc,alpha_SEM);
  boil::oout<<"N_eddy= "<<sem.n_eddy()<<" sig_rep "<<sig_rep<<"\n";
  boil::oout<<"recommendation dt based on sigma= "<<0.2*sig_min/Uc<<"\n";


  /*-------------------+
  |  check if restart  |
  +-------------------*/
  bool restart = false;
  std::fstream input;
  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    restart=true;
  }

  int ts=0;
  if( restart ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
  }
 
  if( restart ) {
    uvw   .load("uvw", ts);
    press .load("press", ts);
#if 1
    /* normal restart */
    wd    .load("wd",   ts);
    // SEM
    // read backup  
    sem.load("sem", ts);  
    // if you change alpha and restart, then
    //sem.initEddies(); // create eddies
    //sem.plot();
    boil::oout<<"N_eddy= "<<sem.n_eddy()<<"\n";
#else
    /* after changeGrid */
    wd    .load("wd",   ts);  // load as an initial value
    sdummy = 0.0;             // reset source term
    di.compute();
    // SEM initialize
    sem.initEddies(); // create eddies
    sem.plot();
    boil::oout<<"N_eddy= "<<sem.n_eddy()<<"\n";
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
    /* compute distance function */
    sdummy = 0.0; // reset source term for di
    di.compute();
    //boil::plot->plot(wd,"wd",0);

    /* velocity */
    Comp m = Comp::u();
    for_vmijk(uvw,m,i,j,k){
      int K = dom.global_K(k)-boil::BW;
      uvw[m][i][j][k] = U_global[K];
    }

    /* SEM: initialize eddies */
    sem.initEddies();
    sem.plot();

    /* plot initial condition */
    boil::plot->plot(uvw,press,mu_t,yplus,"uvw-press-mut-yplus",0);
  }
  input.close();

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;
  int iint2 = int(time.current_time()/tint2) + 1;

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

    /*-----------------+
    |  inlet velocity  |
    +-----------------*/
    alloc3d( &copy_planeVec, 3, uvw.nj()+1, uvw.nk()+1);
    //SEM: convect & recycle eddies
    sem.advance(time.dt());

    //SEM: calculate inlet
    for_m(m){
      for_vmjk(uvw,m,j,k) {
        const real zz = uvw.zc(m,k);
        const real yy = uvw.yc(m,j);
        const int KC = dom.global_K(k)-boil::BW;

        // raw sums s_i(y,z)
        real s_raw[3];
        sem.accumulate(yy, zz, s_raw);  // return: s_raw

        // Final inlet: mean + fluctuation
        if (m==Comp::u()) {
          copy_planeVec[~m][j][k] = U_global[KC] + L11_global[KC]*s_raw[0];
        } else if (m==Comp::v()) {
          copy_planeVec[~m][j][k] =                L22_global[KC]*s_raw[1];
        } else {
          copy_planeVec[~m][j][k] =                L31_global[KC]*s_raw[0]+L33_global[KC]*s_raw[2];
        }
          // p[boil::BW][j][k]=s_raw[0]; // store s_raw[0] to visualize
      }
    }
    // copy_planeVec is computed in all processes but copied to uvw only at inlet
    uvw.bnd_insert(Dir::imin(),copy_planeVec);
    //boil::plot->plot(uvw,"uvw",1);
    //boil::plot->plot(p,"p",1); exit(0);
    //dealloc3d(&copy_planeVec);

    /*---------------------+
    |  reset source terms  |
    +---------------------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /*------------------+
    |  turbulence model |
    +------------------*/
    /* Smagorinsky model */
    // without damping function
    //tm.smagorinsky( & ns, & mu_t, 0.17);
    // with damping function
    //const real dist_max = 1.0e-3; // at z = 1mm, yplus ~ 180. damping ~ 0 if yplus > 40.
    //tm.smagorinsky( & ns, & mu_t, 0.17, &dist_max, &wd, &yplus);  // with damping
    /* WALE model */
    tm.wale( & ns, & mu_t, 0.325);  // 0.325^2 ~ 0.544^2

    /* limit eddy viscosity */
    for_avijk(mu_t,i,j,k) {
      mu_t[i][j][k] = std::min(mu_t[i][j][k],mu_t_max);
    }

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* discretization */
    ns.discretize( &mu_t );
    tm.tau_wall( & ns, wd, & xyz );  // tau_wall must be called after ns.discretize because it accesses A[]
    pr.discretize();
    pr.coarsen();

    /* intermediate velocity */
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1.0e-32));

    /* pressure poisson equation */
    p = 0.0;
    sol.vcycle(ResRat(1.0e-16));

    /* update velocity and pressure */
    ns.project(p);
    press += p;

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(dom.ibody().on(i,j,k)){
        if(pmin>press[i][j][k]) pmin=press[i][j][k];
      }
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k){
      if(dom.ibody().on(i,j,k)){
        press[i][j][k] -= pmin;
      } else {
        press[i][j][k] = 0.0;
      }
    }
    press.bnd_update();
    press.exchange_all();

#if 1
    /* limit velocity */
    for_m(m)
      for_avmijk(uvw,m,i,j,k) {
        uvw[m][i][j][k] = max(-25.0,min(25.0,uvw[m][i][j][k]));
      }
#endif

    /*-------------+
    |  dt control  |
    +-------------*/
    real cflmax = ns.cfl_max();
    boil::oout<<"main:cflmax= "<<time.current_time()<<" "<<cflmax<<"\n";
    /* interface is always included because of outlet */
    time.control_dt(cflmax, cfl_limit, dt_max);

    if (time.dt()<1e-8) {
      boil::oout<<"Too small time step: "<<time.dt()<<"\n";
      exit(0);
    }

    /*--------------+
    |  output data  |
    +--------------*/
    /* tecplot files */
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,press,mu_t,yplus,
                      "uvw-press-mut-yplus",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /* monitoring */
#ifdef OPT_RACK
    #include "monitor_output.cpp"
#endif

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % nint==0) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      wd   .save("wd",   time.current_step());
      sem  .save("sem",  time.current_step());
    }


    if( boil::timer.current_min() > (wmin)
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      wd   .save("wd",   time.current_step());
      sem  .save("sem",  time.current_step());

      uvw  .rm("uvw",   ts);
      press.rm("press", ts);
      wd   .rm("wd", ts);
      sem  .rm("sem", ts);
    }

    if((time.current_step()) % (nint)==0 ) {
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
 
    if( boil::timer.current_min() > (wmin)
      || time.current_step()==time.total_steps()) {
      std::fstream output;
      output << std::setprecision(16);
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      boil::timer.stop();
      boil::timer.report();
      exit(0); 
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}
/*-----------------------------------------------------------------------------+
 '$Id: main-dam.cpp,v 1.12 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
