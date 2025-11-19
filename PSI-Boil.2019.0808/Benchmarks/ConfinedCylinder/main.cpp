#include "Include/psi-boil.h"

#define LEVEL  2 /* 1, 2, 3 */
#define RE   100 /* 20, 100 */

const real H  = 0.41;
const real L  = 2.2;

const real D   = 0.1;
const real MU  = 0.001;
const real RHO = 1.0;

#if LEVEL==1
  const int NY =  64;
#elif LEVEL==2
  const int NY = 128;
#else
  const int NY = 256;
#endif
const int NX1 = NY;
const int NX2 = NY;
const int NZ  = 4;

#if RE==20
  char u_in[] = "4.0*0.3*y*(0.41-y)/0.41^2";
#else
  char u_in[] = "4.0*1.5*y*(0.41-y)/0.41^2";
#endif

#include "boxed_forces.cpp"

Range<int> box_i(-NY/4, NY/4);
Range<int> box_j(-NY/4, NY/4);

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC(AsNodes::yes());

  /*----------------+
  |  immersed body  |
  +----------------*/
  Body cyl("cylinder_bounded.stl");

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gz ( Range<real>(-0.5, 0.5), NZ, Periodic::yes());
  Grid1D gy ( Range<real>(0.0, H), NY, Periodic::no());
  Grid1D gx1( Range<real>(0.0, H), NX1, Periodic::no());

  const real dx = 0.41/(real)NY;

  Grid1D gx2( Range<real>(H, L), 
              Range<real>(dx, 8.0*dx),
              NX2, Periodic::no());
  Grid1D gx( gx1, gx2, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz, &cyl);

  /*----------+ 
  |  ndt, dt  |
  +----------*/
  int ndt = 20000;   
  if(LEVEL==2) ndt *=2;
  if(LEVEL==3) ndt *=4;
  real dt = 0.004;
  if(LEVEL==2) dt /=2;
  if(LEVEL==3) dt /=4;
  if(RE==100) dt /= 5.0;

  Times time(ndt, dt); 

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * smoother = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 
                   u_in, "0.0", "0.0") );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
  }

  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu (MU);
  fluid.rho(RHO);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, smoother, &fluid);
  ns.convection_set(ConvScheme::central());
  //ns.diffusion_set(TimeScheme::backward_euler());
  //ns.convection_set(TimeScheme::forward_euler());
  //ns.discretize();
  //ns.convection_set(ConvScheme::upwind());

  Pressure pr(p, f, uvw, time, smoother, &fluid);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);

  Location loc_0("monitor_0", d, d.I(0.15), d.J(0.27), NZ/2);
  Location loc_1("monitor_1", d, d.I(0.35), NY/2, NZ/2);
  Location loc_2("monitor_2", d, d.I(0.45), NY/2, NZ/2);
  Location loc_3("monitor_3", d, d.I(0.55), NY/2, NZ/2);
  Location loc_4("monitor_4", d, d.I(0.65), NY/2, NZ/2);
  Location loc_5("monitor_5", d, d.I(0.75), NY/2, NZ/2);

  /*------------+
  |  time-loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
/*============================================================================*/
/* I */ const int I = NY/2;
/* I */ const int J = NY/2;
/* I */ const int K =  2;
/* I */
/* I */ const real m_x_0 = innertia_x(uvw, I,J,K, box_i, box_j);
/* I */ const real m_y_0 = innertia_y(uvw, I,J,K, box_i, box_j);
/* I */
/* I */ real fc_x_0 = convective_force_x(uvw, I,J,K, box_i, box_j);
/* I */ real fc_y_0 = convective_force_y(uvw, I,J,K, box_i, box_j);
/*============================================================================*/

    ns.cfl_max();

    ns.new_time_step();
    ns.convection();
    ns.solve(ResRat(1e-2));

/*============================================================================*/
/* I */ const real m_x_1 = innertia_x(uvw, I,J,K, box_i, box_j);
/* I */ const real m_y_1 = innertia_y(uvw, I,J,K, box_i, box_j);
/* I */ const real fi_x_10 = (m_x_1-m_x_0)/time.dt(); /* [kg m / s^2] = [N] */
/* I */ const real fi_y_10 = (m_y_1-m_y_0)/time.dt(); /* [kg m / s^2] = [N] */
/* I */
/* I */ real fv_x_1 = viscous_force_x(uvw, I,J,K, box_i, box_j);
/* I */ real fv_y_1 = viscous_force_y(uvw, I,J,K, box_i, box_j);
/*============================================================================*/

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-2));
    p.exchange();
    ns.project(p);
    pr.update_rhs();

/*============================================================================*/
/* I */ const real m_x_2 = innertia_x(uvw, I,J,K, box_i, box_j);
/* I */ const real m_y_2 = innertia_y(uvw, I,J,K, box_i, box_j);
/* I */
/* I */ const real fi_x_21 = (m_x_2-m_x_1)/time.dt(); /* [kg m / s^2] = [N] */
/* I */ const real fi_y_21 = (m_y_2-m_y_1)/time.dt(); /* [kg m / s^2] = [N] */
/* I */
/* I */ const real pr_x_2 = pressure_force_x(p, I,J,K, box_i, box_j);
/* I */ const real pr_y_2 = pressure_force_y(p, I,J,K, box_i, box_j);
/* I */
/* I */ const real fi_x_20 = (m_x_2-m_x_0)/time.dt(); /* [kg m / s^2] = [N] */
/* I */ const real fi_y_20 = (m_y_2-m_y_0)/time.dt(); /* [kg m / s^2] = [N] */
/*============================================================================*/
     
//  OMS(---x---);
//  OPR(fc_x_0+fv_x_1        - fi_x_10);
    const real fx = fc_x_0+fv_x_1+pr_x_2 - fi_x_20;
    OPR( fx );
//  OMS(---y---);
//  OPR(fc_y_0+fv_y_1        - fi_y_10);
    const real fy = fc_y_0+fv_y_1+pr_y_2 - fi_y_20;
    OPR( fy );

    #if RE==20
      const real DIV = 0.004; 
    #else
      const real DIV = 0.1; 
    #endif

    boil::oout << "c_D_new = " << 2.0 * fx / DIV << boil::endl;
    boil::oout << "c_L_new = " << 2.0 * fy / DIV << boil::endl;

    /*-----------------------------+
    |  print monitoring locations  |
    +-----------------------------*/
    loc_0.print(p);
    loc_1.print(uvw, Comp::v());
    loc_2.print(uvw, Comp::v());
    loc_3.print(uvw, Comp::v());
    loc_4.print(uvw, Comp::v());
    loc_5.print(uvw, Comp::v());

    /*------------+
    |  save/plot  |
    +------------*/
    int period = 1000;
    if(LEVEL==2) period*=2;   
    if(LEVEL==3) period*=4;   

    if(time.current_step() % period == 0) {
      std::stringstream u_name_bck;
      std::stringstream p_name_bck;

      if(RE==20) {
        u_name_bck << "uvw_L" << LEVEL << "_R020";
        p_name_bck << "p_L"   << LEVEL << "_R020";
      } else {
        u_name_bck << "uvw_L" << LEVEL << "_R100";
        p_name_bck << "p_L"   << LEVEL << "_R100";
      }

      uvw.save(u_name_bck.str().c_str(), time.current_step());
      p.  save(p_name_bck.str().c_str(), time.current_step());
    }

    if(time.current_step() % period == 0) {
      std::stringstream up_name_dat;
      std::stringstream ns_name_bck;
      std::stringstream pr_name_bck;

      up_name_dat << "uvw_L" << LEVEL;
      ns_name_bck  << "ns_L" << LEVEL;
      pr_name_bck  << "pr_L" << LEVEL;
      if(RE==20) {up_name_dat << "_R020"; 
                  ns_name_bck << "_R020"; 
                  pr_name_bck << "_R020";}
      else       {up_name_dat << "_R100"; 
                  ns_name_bck << "_R100"; 
                  pr_name_bck << "_R100";}
      boil::plot->plot(uvw, p, up_name_dat.str().c_str(), time.current_step());
      ns.save(ns_name_bck.str().c_str(), time.current_step());
      pr.save(pr_name_bck.str().c_str(), time.current_step());
    }

    /*---------+
    |  exit ?  |
    +---------*/
    std::ifstream infile;
    infile.open ("stop.now", std::ifstream::in);
    if(infile.good()) exit(0);
    infile.close();
  }
  boil::plot->plot(uvw,  p, "uvw,p",  time.current_step()-1);
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
