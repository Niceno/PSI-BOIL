#include "Include/psi-boil.h"

/* domain parameters */
const real LX = 8.0e-3;
//const int  NX = 8;
const int  NX = 64;
const int  NY = 4;
const real LY = LX/real(NX)*real(NY);

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
  Grid1D gx( Range<real>(0,LX), 
             NX, Periodic::no());

  Grid1D gy( Range<real>(0,LY),  
             NY, Periodic::no(),
             BndGrid::symmetry());  // Yohei

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gy);

  /*-------+
  |  time  |
  +-------*/
  const int  ndt = 1000;
  const real dt = 0.001; 
  Times time(ndt, dt);
  const int nint = 100; // interval of field-data output
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // pressure
  Scalar h  (d), g  (d); // enthalpy    
  Scalar t  (d);         // temperature
  Scalar rho_old (d);    // density before update

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  f = p.shape(); // Yohei
  g = p.shape(); // Yohei
  t = p.shape();
  rho_old = p.shape();

  h.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 133780 ) ); // @ -30
  h.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), 284300 ) ); // @ +30.05
  h.bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
  h.bc().add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
  h.bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
  h.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  LookUpTable CO2("co2.txt");

  Matter fluid(d);
  fluid.mu (0.078078);
  fluid.rho(700.72);
  fluid.cp (5255.8);

  fluid.variable(Set::rho()); /* fluid is variable and a function of t */
  fluid.variable(Set::mu());  /* fluid is variable and a function of t */

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr  (p,   f,   uvw, time, solver, &fluid);
  Momentum ns  (uvw, xyz,      time, solver, &fluid);
  Enthalpy enth(h,   g,   uvw, time, solver, &fluid);
  enth.diffusion_set(TimeScheme::backward_euler());

  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  h = 284300;  // applied to dummy cells as well; exchange is unnecessary
  t.look_up(h, CO2, Column(2), Column(0));  // applied to dummy cells; exchange is unnecessary
  fluid.look_up(Set::rho(), t, CO2, Column(0), Column(3)); // applied to dummy cells
  fluid.look_up(Set::mu(),  t, CO2, Column(0), Column(6)); // applied to dummy cells
  boil::plot->plot(uvw,p,h,t,*fluid.rho(),"uvw-p-h-t-rho",0);

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /* store density to rho_old */
#if 1
    /* loop including dummy cells */
    for_avijk(rho_old,i,j,k){
      rho_old[i][j][k] = fluid.rho(i,j,k);
    }
#else
    for_vijk(rho_old,i,j,k){
      rho_old[i][j][k] = fluid.rho(i,j,k);
    }
    rho_old.exchange();
#endif

    /*---------------+
    |  enthalpy eq.  |
    +---------------*/
    enth.discretize();
    enth.new_time_step();
    enth.solve(ResRat(1e-8), "enth");
    /* update temperature, rho and mu */
    t.look_up(h, CO2, Column(2), Column(0)); // update temperature
    fluid.look_up(Set::rho(), t, CO2, Column(0), Column(3));// update rho
    fluid.look_up(Set::mu(),  t, CO2, Column(0), Column(6));// update mu

    /*---------------+
    |  momentum eq.  |
    +---------------*/
    /* set rhs for pressure poisson equation */
    for_vijk(f,i,j,k){
      f[i][j][k]=-1.0/time.dt()*f.dV(i,j,k)/fluid.rho(i,j,k)
                 *(fluid.rho(i,j,k)-rho_old[i][j][k])/time.dt();
    }
    f.exchange();
    ns.vol_phase_change(&f); // set velocity at outlet to conserve mass

    ns.discretize();
    pr.discretize();
    pr.coarsen();
    ns.new_time_step();
    ns.solve(ResRat(0.0001));
    p = 0.0;
    multigrid.vcycle(ResRat(1e-3));
    p.exchange();
    ns.project(p); // update velocity
    ns.cfl_max();  // output cfl-max

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw,p,h,t,*fluid.rho(),
            "uvw-p-h-t-rho",time.current_step());
    }
    /* warning: next part works only single process (non-paralell) */
    if(time.current_step() == 1) {
      std::fstream output;
      output << std::setprecision(16);
      output.open("field-1st-step.txt", std::ios::out);
      int j=t.sj();
      int k=t.sk();
      for_vi(t,i){
        output<<i<<" "<<t[i][j][k]<<" "<<fluid.rho(i,j,k)<<" "<<rho_old[i][j][k]<<" "
              <<f[i][j][k]<<" "<<uvw[Comp::u()][i][j][k]<<"\n";
      }
      output.close();
    }

    boil::oout<<"main: min&max: "<<time.current_time()<<" "<<p.min()<<" "<<p.max()<<" "
              <<h.min()<<" "<<h.max()<<" "<<t.min()<<" "<<t.max()<<"\n";
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();
}

