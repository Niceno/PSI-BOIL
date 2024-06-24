#include "Include/psi-boil.h"

/* domain parameters */
const real LX = 1.0;
const real LY = 0.125;
const int  NX = 64;
const int  NY = 4;

const real t_init=50.0;
const real dT = 0.001521;

/******************************************************************************/
int main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no());
  Grid1D gy( Range<real>(0,LY), NY, Periodic::yes()); 

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);         // vel, RHS of momentum
  Scalar p  (d), press(d), f(d); // dp, pressure, RHS of pressure
  Scalar t  (d), q  (d);         // temperature, RHS of temperature
  Scalar rho_old (d);            // density before update

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  f = p.shape();
  q = p.shape();
  rho_old = p.shape();
  press = p.shape();

  t.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), t_init+0.5*dT ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), t_init-0.5*dT ) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  LookUpTable air("air3.txt");
  boil::oout<<"success to open LookUpTable.\n";

  const real rho_init = 1.0816667;    // at 50 deg.
  const real mu_init = 0.0000197;     // at 50 deg.
  const real lambda_init = 0.0277667; // at 50 deg.
  const real cp_init = 1092.672659;   // at 50 deg.
  Matter fluid(d);
  fluid.rho(rho_init);
  fluid.mu(mu_init);
  fluid.lambda(lambda_init);
  fluid.cp(cp_init); // J/m3.K
  const real beta = 0.0031050;  // at 50 deg. used only for Ra calculation
  const real gravity = 9.8;
  const real alpha = fluid.lambda()->value()/fluid.cp()->value();
  const real nu = fluid.mu()->value()/fluid.rho()->value();
  boil::oout<<"alpha, beta, nu= "<<alpha<<" "<<beta<<" "<<nu<<"\n";
  boil::oout<<"Pr= "<<nu/alpha<<"\n";
  boil::oout<<"Ra= "<<fluid.rho()->value()*gravity*beta*dT
                     /(fluid.mu()->value()*alpha)<<"\n";

  fluid.variable(Set::rho()); /* fluid is variable and a function of t */
  fluid.variable(Set::mu());  /* fluid is variable and a function of t */
  fluid.variable(Set::cp());  /* fluid is variable and a function of t */
  fluid.variable(Set::lambda());  /* fluid is variable and a function of t */

  /*-------+
  |  time  |
  +-------*/
  const int  ndt = 4000;
  const real dt = 0.000075/alpha;
  const real cfl_limit = 0.25;
  Times time(ndt, dt);
  const int nint = 1000; // interval step of field-data output
  const real tint = 1000.0; // interval time of field-data output

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::di());
  Pressure pr  (p,   f,   uvw, time, solver, &fluid);
  Momentum ns  (uvw, xyz,      time, solver, &fluid);
  Enthalpy enth(t,   q,   uvw, time, solver, &fluid);
  enth.diffusion_set(TimeScheme::backward_euler());

  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  t = t_init;  // applied to also dummy cells as well; exchange is unnecessary
  fluid.look_up(Set::rho(),    t, air, Column(0), Column(3)); // applied to dummy cells
  fluid.look_up(Set::cp(),     t, air, Column(0), Column(4)); // applied to dummy cells
  fluid.look_up(Set::lambda(), t, air, Column(0), Column(5)); // applied to dummy cells
  fluid.look_up(Set::mu(),     t, air, Column(0), Column(6)); // applied to dummy cells
  /* hydro-static pressure */
  for_vijk(press,i,j,k){
    press[i][j][k]= -gravity * rho_init * press.zc(k);
  }
  press.exchange();
  /* output initial condition */
  boil::plot->plot(uvw,press,t,*fluid.rho(),*fluid.lambda(),*fluid.cp(),*fluid.mu(),
            "uvw-press-t-rho-lambda-cp-mu",0);
  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /* store density to rho_old */
    /* loop including dummy cells */
    for_avijk(rho_old,i,j,k){
      rho_old[i][j][k] = fluid.rho(i,j,k);
    }

    /*---------------+
    |  enthalpy eq.  |
    +---------------*/
    enth.discretize();
    enth.new_time_step();
    enth.solve(ResRat(1e-8), "enth");

    /* update rho, cp, lambda and mu */
    fluid.look_up(Set::rho(),    t, air, Column(0), Column(3)); // applied to dummy cells
    fluid.look_up(Set::cp(),     t, air, Column(0), Column(4)); // applied to dummy cells
    fluid.look_up(Set::lambda(), t, air, Column(0), Column(5)); // applied to dummy cells
    fluid.look_up(Set::mu(),     t, air, Column(0), Column(6)); // applied to dummy cells

    /*---------------+
    |  momentum eq.  |
    +---------------*/
    /* gravity force */
    const Comp m = Comp::w();
    for_avmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = -gravity * xyz.dV(m,i,j,k) * fluid.rho(m,i,j,k);

#if 1
    /* set rhs for pressure poisson equation related to density change*/
    for_vijk(f,i,j,k){
      f[i][j][k]=-1.0/time.dt()*f.dV(i,j,k)/fluid.rho(i,j,k)
                 *(fluid.rho(i,j,k)-rho_old[i][j][k])/time.dt();
    }
    f.exchange();

    /* try to satisfy volume conservation */
    real sum_f = 0.0;
    int sum_cell = 0;
    for_vijk(f,i,j,k){
        sum_f += f[i][j][k];
        sum_cell++;
    }
    boil::cart.sum_real(&sum_f);
    boil::cart.sum_int(&sum_cell);
    real ave_f = sum_f/real(sum_cell);
    for_vijk(f,i,j,k){
      f[i][j][k] -= ave_f;
    }
    f.exchange();

    /* set velocity at outlet to conserve mass */
    ns.vol_phase_change(&f);
#endif

    ns.discretize();
    pr.discretize();
    pr.coarsen();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(0.0001));
    p = 0.0;
    multigrid.vcycle(ResRat(1e-3));
    p.exchange();
    ns.project(p); // update velocity
    for_vijk(press,i,j,k){
      press[i][j][k] += p[i][j][k];
    }
    press.exchange();

    /*-------------+
    |  dt control  |
    +-------------*/
    real cflmax=ns.cfl_max();
#if 0
    if(cflmax>1.0){
      boil::oout<<"main: Too large CFL!!!";
      exit(0);
    }
#else
    time.control_dt(cflmax, cfl_limit, dt);
#endif

    /*-------------+
    |  field data  |
    +-------------*/
#if 0
    if((time.current_time()) / (tint) >= real(iint) ) {
      boil::plot->plot(uvw,press,t,*fluid.rho(),*fluid.lambda(),*fluid.cp(),*fluid.mu(),
            "uvw-press-t-rho-lambda-cp-mu",time.current_step());
      iint++;
    }
#endif
#if 1
    if(time.current_step()%nint ==0 ) {
      boil::plot->plot(uvw,press,t,*fluid.rho(),*fluid.lambda(),*fluid.cp(),*fluid.mu(),
            "uvw-press-t-rho-lambda-cp-mu",time.current_step());
    }
#endif

    boil::oout<<"main: min&max: "<<time.current_time()<<" "<<t.min()<<" "<<t.max()<<"\n";
  }

  /* output field data: (field at the end) - (initial value) */
  for_vijk(press,i,j,k){
    t[i][j][k] -= t_init;
    rho_old[i][j][k]=fluid.rho(i,j,k)-rho_init;
    p[i][j][k]=fluid.lambda(i,j,k)-lambda_init;
    q[i][j][k]=fluid.cp(i,j,k)-cp_init;
    f[i][j][k]=fluid.mu(i,j,k)-mu_init;
  }
  boil::plot->plot(uvw,press,t,rho_old,p,q,f,
            "uvw-press-dt-drho-dlambda-dcp-dmu",time.current_step()-1);

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();
}
