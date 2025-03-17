/*----------------------------------------------------------------------------+
|  Rising single CO2 bubble gas in stagnant water                             |
|  using reference frame                                                      |
|  for publication: Direct numerical simulation of a single CO2 bubble rising |
|                   and dissolving in quiescent water                         |
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"
#include <vector>
#include <iomanip>  // For precision control
#include <cmath>
//#define STOP_MASS_TRANSFER  // activate: no mass transfer
//#define USE_SZPLT  // activate: output szplt, otherwise ascii tecplot dat
using namespace std;

// Case a
#if 0
const real dia_target = 2.176e-4;
const real surface_tension_coeff = 4.5939e-4;
#ifdef STOP_MASS_TRANSFER
  const real dia_init = dia_target;
#else
  const real dia_init = dia_target * 1.6;  // initial bubble diameter
#endif
#endif

// Case b
#if 0
const real dia_target = 0.00064781;  // b
const real surface_tension_coeff = 1.0179E-4; // b
#ifdef STOP_MASS_TRANSFER
  const real dia_init = dia_target;
#else
  const real dia_init = dia_target * 1.3;  // initial bubble diameter
#endif
#endif

// Case c
#if 0
const real dia_target = 0.00093034;  // c
const real surface_tension_coeff = 2.6872E-3; // c
#ifdef STOP_MASS_TRANSFER
  const real dia_init = dia_target;
#else
  const real dia_init = dia_target * 1.2;  // initial bubble diameter
#endif
#endif

// Case CO2
#if 1
const real dia_target = 0.009;  // [m]
const real surface_tension_coeff = 0.07; // liquid water - CO2 gas
const real dia_init = dia_target;
#endif

const real rad_target = 0.5 *dia_target;
const real rad_init   = 0.5 *dia_init;


const int NX= 2*64;
const int NZ= NX*1;
const real LX = 10*rad_target; // width = 5*db m
const real LZ = LX*(NZ/NX);

const real gravity = -9.8; // [m/s2]
const real c_surf = 1.0; // [kg/m3] -> non-dimensionalized value

const real zcent = LZ - 2.0*rad_target - rad_init;  // z of initial bubbles [m]
const real z_bubble_top_target = LZ - 2.0*rad_target;  // target z of bubble top [m]

vector<real> solveLeastSquares(const vector<vector<real>>& A, const vector<real>& b);
double computeVelocityLeastSquares(const vector<double>& timeLS, const vector<double>& heightLS, int poly_order = 2);

/******************************************************************************/
int main(int argc, char * argv[]) {

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
#ifndef USE_SZPLT
  boil::plot = new PlotTEC();
#else
  boil::plot = new PlotTECMPI();
#endif


  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), NX, Periodic::no());
  Grid1D gz( Range<real>(0.0, LZ), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter gas(d), water(d);
  gas.rho   (1.977);
  gas.mu    (1.48e-5);
  gas.lambda(1.88e-9);
  water.rho   (997.0);
  water.mu    (1.0e-3); //dynamic viscosity [Pa.s]
  water.lambda(1.88e-9); //D [m2/s]

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d), press(d,"press"); // pressure
  Scalar c  (d,"vfl"), g  (d), kappa(d); // concentration
  Scalar tpr(d,"eps"), q  (d); // 
  //Scalar idFlood(d,"idFlood");
  Scalar mdot(d,"mdot");  // phase change rate [kg/m3s]
                          // negative for dissolution of CO2 (gas) into liquid water
  vector<real> timeLS(10,0.0), height(10,0.0);  // for bubble velocity calculation
                  
  real w_inlet = 0.0;  // inlet velocity [m/s]
  // camera = origin of coordinate system
  real w_camera = -w_inlet; // camera velocity in Z [m];
  real z_camera = 0.0; // camera position in Z [m]

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::inlet(), 0.0, 0.0, w_inlet ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  press = p.shape();
  f = p.shape();
  g = p.shape();
  kappa = p.shape();
  //idFlood = p.shape();
  mdot = p.shape();
  
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 1.0 ) );

  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::outlet() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );

  Matter mixed(water, gas, & c);
  mixed.sigma(surface_tension_coeff);

  boil::oout<<"Schmidt number = "<<water.mu()->value()/water.rho()->value()/water.lambda()->value()<<"\n";
  const real delta_rho = water.rho()->value() -gas.rho()->value();
  boil::oout<<"Eotvos number will decrease due to disolution!\n";
  boil::oout<<"Eotvos number = "<<fabs(gravity)*delta_rho*pow(dia_target,2.0)/mixed.sigma()->value()<<"\n";
  boil::oout<<"Morton number = "<<fabs(gravity)*pow(water.mu()->value(),4.0)*delta_rho/(pow(water.rho()->value(),2.0)*pow(mixed.sigma()->value(),3.0))<<"\n";

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 5.0 * pow(gas.rho()->value()*pow(dxmin,3.0)
                        /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 10000000;
  Times time(ndt, dt);
  const real tint = 1.0e-3;  // time interval for Tecplot output
  const int nint= 10000;      //set the backup time interval
  const real cfl_limit=0.25; // max of CFL

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  /* linear solver */
  Krylov * solver = new CG(d, Prec::ic2());

  /* Navier-Stokes */
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.set_min_iteration(3);
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.max_cycles(10);
  multigrid.min_cycles(4);

  /* VOF */
  VOF conc (c,   g, kappa, uvw, time, solver);

  /* enthalpy equation */
  EnthalpyFDhighPr enthFD(tpr, q, c, uvw, time, solver, & mixed ,c_surf);
  //EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, & mixed ,c_surf);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  enthFD.set_conv_divu_subtract(false);  // Conservative form for spieces equation nabla(cu)

  /* phase change */
  const real latent = 1.0;
  PhaseChange pc(mdot, tpr, q, c, g, f, c, uvw, time, & mixed, latent, c_surf);
  pc.set_rhol(1.0e+12);
  //pc.set_turbP(1.0);  // activate when turbulent eddy is intriduced
 
  /* floodfill */
  //Floodfill flood(c, idFlood, &uvw, time);
  //flood.set_out_freq(1);
  //boil::oout<<"main:Floodfill:output frequency= "<<flood.get_out_freq()<<"\n";
        
  /* PIDcontrol */
  PIDcontrol pid_control(1e+6,0.0,1e+4);  // coeffs for Proportional (=1e+6), Integral and Derivative

  /*-------------------------------+
  |  initial condition or restart  |
  +-------------------------------*/
  int ts=0;
  bool restart = false;

  std::fstream input;
  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    restart=true;
  }

  if( restart ) {
    /*----------+
    |  restart  |
    +----------*/
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    input >> w_inlet;
    input >> z_camera;
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);

    uvw.load     ("uvw",   ts);
    press.load   ("press", ts);
    c.load       ("c",     ts);
    tpr.load     ("tpr",   ts);
    pid_control.load("pidcont",ts);  // load information for PID controller
    fstream input_velb;
    input_velb.open("velb.bck",std::ios::in);
    for(int i = 0; i<10; i++) {
      input_velb>>timeLS[i]>>height[i];
      boil::oout<<"main:velb "<<i<<" "<<timeLS[i]<<" "<<height[i]<<"\n";
    }
    input_velb.close();

    // update boundary condition
    for_m(m) {
      uvw.bc(m).modify( BndCnd( Dir::imin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::imax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::jmin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::jmax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::kmax(), BndType::inlet(), 0.0, 0.0, w_inlet ) );
    }
    conc.init();  // This is mandatory to detect the interface position at the beginning of time loop
    conc.front_minmax();

  } else {
    /*--------------------+
    |  initial condition  |
    +--------------------*/
    Comp m = Comp::w();
    for_vmijk(uvw,m,i,j,k){
      uvw[m][i][j][k] = w_inlet;
    }

    for_vijk(c,i,j,k)
      c[i][j][k] = 1.0;
  
    for_vijk(c,i,j,k) {
      real dist=sqrt(pow(c.xc(i),2.0)+pow(c.yc(j),2.0)+pow((c.zc(k)-zcent),2.0));
      if (dist<rad_init*0.75) {
        c[i][j][k]=0.0;
      } else if(dist<rad_init*1.25) {
        int mm=8;
        real x0=d.xn(i);
        real y0=d.yn(j);
        real z0=d.zn(k);
        real ddx=d.dxc(i)/real(mm);
        real ddy=d.dyc(j)/real(mm);
        real ddz=d.dzc(k)/real(mm);
        int itmp=0;
        for (int ii=0; ii<mm; ii++){
          for (int jj=0; jj<mm; jj++){
            for (int kk=0; kk<mm; kk++){
              real xxc=x0+0.5*ddx+real(ii)*ddx;
              real yyc=y0+0.5*ddy+real(jj)*ddy;
              real zzc=z0+0.5*ddz+real(kk)*ddz;
              real dist=sqrt(pow(xxc,2.0)+pow(yyc,2.0)+pow(zzc-zcent,2.0));
              if (dist>rad_init){
                itmp=itmp+1;
              }
            }
          }
        }
        c[i][j][k]=real(itmp)/real(mm*mm*mm);
      }
    }
    c.bnd_update();
    c.exchange_all();

    conc.init();  // This is mandatory  2024.07.05
    conc.front_minmax();

    for_vijk(c,i,j,k) {
      if (c[i][j][k]>0.5) {
        tpr[i][j][k] = 0.0;
      } else {
        tpr[i][j][k] = c_surf;
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();

#ifndef USE_SZPLT
    boil::plot->plot(uvw,c,press,tpr,mdot,"uvw-c-press-tpr-mdot",0);
#else
    boil::plot->plot("uvw-c-press-tpr-mdot",0,&time,&uvw,&c,&press,&tpr,&mdot);
#endif

  }
  input.close();

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  Time loop  |
  +------------*/

  for(time.start(); time.end(); time.increase()) {

    /* reset body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* change w_inlet */
    real z_bubble_top_current = conc.get_zmaxft(); // z of bubble top [m]
#if 0
    real coef_z = 1.0e+8;  // coefficient of camera work [1/s2], how quickly follow the bubble
    // acceleration of camera
    real acc_camera = coef_z * (z_bubble_top_current - z_bubble_top_target); // [m/s2]
    real acc = -acc_camera;     // acceleration of inlet velocity [m/s2]
#endif
#if 1
    real setpoint = z_bubble_top_target;
    real measured = z_bubble_top_current;
    // - acceleration of coordinate system
    real acc = pid_control.calculate(setpoint,measured,time.current_time());
    real acc_camera = -acc;     // acceleration of coordinate system (= camera) [m/s2]
#endif

    // bubble volume center height
    real volb=0.0;
    real volb_height=0.0;
    for_vijk(c,i,j,k){
      volb += (1.0-c[i][j][k]) * c.dV(i,j,k); // [m3]
      volb_height += (1.0-c[i][j][k]) * c.dV(i,j,k) * c.zc(k); // [m4]
    }
    boil::cart.sum_real(&volb);
    boil::cart.sum_real(&volb_height);
    real z_bubble_cent = volb_height/volb;  //[m]

    // bubble velocity
    for(int i = 1; i<=9; i++){
      timeLS[i-1]=timeLS[i];
      height[i-1]=height[i];
    }
    timeLS[9]=time.current_time();
    height[9]=z_camera+z_bubble_cent;
    //for(int i = 0; i<=9; i++){
    //  boil::oout<<"main:velocity,i= "<<i<<" "<<timeLS[i]<<" "<<height[i]<<"\n";
    //}
    // bubble velocity computed with least square using 10 past time and position
    real w_bubble = computeVelocityLeastSquares(timeLS, height);

    w_inlet += acc * time.dt();  // inlet velocity [m/s]
    w_camera = -w_inlet;           // camera velocity [m/s]
    z_camera += w_camera * time.dt();  // camera position [m]
    boil::oout<<"main:w_inlet:time "<<time.current_time()<<" w_inlet "<<w_inlet
              <<" accl_camera "<<acc_camera<<" z_camera "<<z_camera
              <<" z_bubble_top_earth "<<z_camera+z_bubble_top_current
              <<" z_bubble_center_earth "<<z_camera+z_bubble_cent
              <<" w_bubble "<<w_bubble<<"\n";
    // update boundary condition
    for_m(m) {
      uvw.bc(m).modify( BndCnd( Dir::imin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::imax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::jmin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::jmax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::kmax(), BndType::inlet(), 0.0, 0.0, w_inlet ) );
    }

    /* phase change */
    pc.update();
#ifdef STOP_MASS_TRANSFER
    mdot = 0.0;  f = 0.0;  g = 0.0;  q = 0.0;
#else
    ns.vol_phase_change(&f);
#endif

    /* body force */
    // gravity and acceleration of reference frame (camera)
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = (gravity+acc)*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }
    xyz.exchange();

    // surface tension
    conc.tension(&xyz, mixed);

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(0.001));

    p = 0.0;
    //multigrid.max_cycles(10);
    multigrid.vcycle(ResRat(1e-3));
    ns.project(p);
    press += p;    // p_n+1 = p_n + dp
    press.exchange();

    /* solve transport equation for volume fraction of liquid */
    conc.new_time_step();
    conc.advance();
    conc.totalvol();
    conc.front_minmax();

    /* solve species transport equation */
    enthFD.discretize();     // mu_t
    enthFD.new_time_step();  // mu_t
    enthFD.solve(ResRat(1e-16),"enthFD");

    // crude operation Cg = const at c_surf
    for_vijk(tpr,i,j,k) {
      if (c[i][j][k] < 0.5) tpr[i][j][k] = c_surf;
    }
    tpr.bnd_update();
    tpr.exchange();

    c.bnd_update();
    c.exchange_all();
    
    /* floodfill */
    //flood.identify_regions();

    /* dt control */
    time.control_dt(ns.cfl_max(), cfl_limit, dt);

    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);

      // subtract w_inlet just for visualization, w is on the earth fixed coordinate system
      m = Comp::w();
      for_avmijk(uvw,m,i,j,k) { uvw[m][i][j][k] += w_camera; }

#ifndef USE_SZPLT
      boil::plot->plot(uvw,c,press,tpr,mdot,"uvw-c-press-tpr-mdot",iint);
#else
      boil::plot->plot("uvw-c-press-tpr-mdot",iint,&time,&uvw,&c,&press,&tpr,&mdot);
#endif

      // return uvw to the velocity reference to the moving coordinate system
      for_avmijk(uvw,m,i,j,k) { uvw[m][i][j][k] -= w_camera; }

      iint = int(time.current_time()/tint) + 1;
    }

    /* output to console */
    // equivalent bubble diameter
    real dia = pow((6.0*conc.get_vol1()/boil::pi),1.0/3.0);

    // bubble area
    real total_area = 0.0;
    for_vijk(c,i,j,k){
      total_area += conc.adens[i][j][k] * c.dV(i,j,k); // [m2/m3]*[m3]
    }
    boil::cart.sum_real(&total_area);

    real total_flux_CL = 0.0;
    real z_examine = 1.0*rad_target;
    for_vijk(tpr,i,j,k) {
      if (tpr.zn(k)<=z_examine && z_examine<tpr.zn(k+1)) {
        total_flux_CL += tpr[i][j][k]*tpr.dSz(i,j,k)*uvw[Comp::w()][i][j][k];
        // [kg/m3]*[m2]*[m/s]=[kg/s]
      }
    }
    boil::cart.sum_real(&total_flux_CL);

    real eo = fabs(gravity)*delta_rho*pow(dia,2.0)/mixed.sigma()->value();
    real re = (w_bubble)*dia/(water.mu()->value()/water.rho()->value());
    real sherwood1 = -pc.get_smdot_neg()/total_area*dia/water.lambda()->value();
    real sherwood2 = total_flux_CL/total_area*dia/water.lambda()->value();
    //real sherwood1 = -pc.get_smdot_neg()/total_area*(2.0*rad_target)/water.lambda()->value();
    //real sherwood2 = total_flux_CL/total_area*(2.0*rad_target)/water.lambda()->value();

    boil::oout<<"mass_transfer_rate:time= "<<time.current_time()<<" dt[s] "<<time.dt()
        <<" diameter[m] "<<dia
        <<" total_area[m2] "<<total_area
        <<" smdot_neg[kg/s] "<<pc.get_smdot_neg()
        <<" total_flux_CL[kg/s] "<<total_flux_CL
        <<" eo "<<eo
        <<" re "<<re
        <<" Sherwood1 "<<sherwood1
        <<" Sherwood2 "<<sherwood2
        <<"\n";

    if( time.current_step() % nint == 0 ) {
      uvw.save     ("uvw",  time.current_step());//save the data as the binary file '.bck'
      press.save   ("press",time.current_step());
      c.save       ("c",    time.current_step());
      tpr.save       ("tpr",    time.current_step());
      pid_control.save("pidcont",time.current_step());
      //pathline.save("pathline",time.current_step());
      if( boil::cart.iam()==0) {
        std::fstream output;
        output << std::setprecision(16);
	std::string name = name_file("time", ".txt", time.current_step());
        output.open(name, std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output << w_inlet << boil::endl;
        output << z_camera << boil::endl;
        output.close();
        name = name_file("velb", ".bck", time.current_step());
        output.open(name, std::ios::out);
        for(int i = 0; i<10; i++) {
          output<<timeLS[i]<<" "<<height[i]<<"\n";
        }
        output.close();
      }
    }
    if(boil::timer.current_min() > wmin
      || time.current_step()==time.total_steps()) {
      uvw.save     ("uvw",  time.current_step());//save the data as the binary file '.bck'
      press.save   ("press",time.current_step());
      c.save       ("c",    time.current_step());
      tpr.save     ("tpr",  time.current_step());
      pid_control.save("pidcont",time.current_step());
      //pathline.save("pathline",time.current_step());
      if( boil::cart.iam()==0) {
        std::fstream output;
        output << std::setprecision(16);
        output.open("time.txt", std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output << w_inlet << boil::endl;
        output << z_camera << boil::endl;
        output.close();
        output.open("velb.bck", std::ios::out);
        for(int i = 0; i<10; i++) {
          output<<timeLS[i]<<" "<<height[i]<<"\n";
        }
        output.close();
      }
      uvw.rm  ("uvw",ts);  // rm bck files previously used for the restart
      press.rm("press",ts);
      c.rm    ("c",ts);
      tpr.rm  ("tpr",ts);
      pid_control.rm("pidcont",ts);
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }

  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}

/******************************************************************************/
// Function to solve (A^T A) x = A^T b using Gaussian elimination
vector<double> solveLeastSquares(const vector<vector<double>>& A, const vector<double>& b) {
    int m = A.size();         // Number of data points (rows)
    int n = A[0].size();      // Number of coefficients (columns)

    // Compute A^T * A (n × n matrix)
    vector<vector<double>> ATA(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < m; k++) {
                ATA[i][j] += A[k][i] * A[k][j];
            }
        }
    }

    // Compute A^T * b (n × 1 vector)
    vector<double> ATb(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < m; k++) {
            ATb[i] += A[k][i] * b[k];
        }
    }

    // Solve (A^T A) x = (A^T b) using Gaussian elimination
    vector<double> x(n, 0.0);
    for (int i = 0; i < n; i++) {
        // Pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(ATA[k][i]) > fabs(ATA[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(ATA[i], ATA[maxRow]);
        swap(ATb[i], ATb[maxRow]);

        // Make diagonal 1
        if (fabs(ATA[i][i]) < 1e-120) {
          boil::oout<<"Singular matrix: Cannot solve. ATA[i][i]="
                    <<i<<" "<<ATA[i][i]<<"\n";
          //throw runtime_error("Singular matrix: Cannot solve.");
          return vector<double>(n, 0.0);  // Return zero vector instead of throwing an error

        }

        for (int k = i + 1; k < n; k++) {
            double factor = ATA[k][i] / ATA[i][i];
            for (int j = i; j < n; j++) {
                ATA[k][j] -= factor * ATA[i][j];
            }
            ATb[k] -= factor * ATb[i];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        x[i] = ATb[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= ATA[i][j] * x[j];
        }
        x[i] /= ATA[i][i];
    }

    return x;
}
/******************************************************************************/
// Function to compute velocity using least squares polynomial fit
double computeVelocityLeastSquares(const vector<double>& timeLS, const vector<double>& height, int poly_order) {
    int n = timeLS.size();
    if (n == 1) return 0.0;  // One data point -> velocity = 0
    if (n == 2) return (height[1] - height[0]) / (timeLS[1] - timeLS[0]);  // Two points -> simple velocity

    // Limit to at most 10 points
    int max_points = 10;
    vector<double> timeLS_trimmed, height_trimmed;
    if (n > max_points) {
        timeLS_trimmed.assign(timeLS.end() - max_points, timeLS.end());
        height_trimmed.assign(height.end() - max_points, height.end());
        n = max_points;
    } else {
        timeLS_trimmed = timeLS;
        height_trimmed = height;
    }

    // Construct Vandermonde matrix A and vector b
    vector<vector<double>> A(n, vector<double>(poly_order + 1, 0.0));
    vector<double> b(n, 0.0);
    for (int i = 0; i < n; i++) {
        double t_power = 1.0;
        for (int j = 0; j <= poly_order; j++) {
            A[i][j] = t_power;
            t_power *= timeLS_trimmed[i];
        }
        b[i] = height_trimmed[i];
    }

    try {
        // Solve for polynomial coefficients using least squares
        vector<double> coeffs = solveLeastSquares(A, b);
        if (coeffs.size() != poly_order + 1) {
            throw runtime_error("Least squares returned incorrect size.");
        }

        // Extract coefficients
        double a0 = coeffs[0];  // Constant term
        double a1 = coeffs[1];  // Linear term
        double a2 = coeffs[2];  // Quadratic term

        // Compute derivative coefficients
        double velocity;
        double t_last = timeLS_trimmed.back();
        velocity = 2 * a2 * t_last + a1;

        return velocity;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
        return 0.0;
    }
}


