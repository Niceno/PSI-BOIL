/*----------------------------------------------------------------------------+
|  Rising single bubble of artificial gas in stagnant water                   |
|  using reference frame                                                      |
|  materiap property etc: K. Hayashi, A. Tomiyama, Interface Tracking         |
|  Simulation of Mass Transfer from a Dissolving Bubble,                      |
|  Journal of Computational Multiphase Flows, 3 (2011) 247-261.               |
|  10.1260/1757-482x.3.4.247.                                                 |
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"

const int NX= 128;   // 2048 is the maximum case
const int NZ= NX*1;
const real radius = 0.5 *0.009;  // diameter = 0.2177 mm
//const real radius_a = radius;
//const real radius_b = radius;
const real radius_a = 5.669644725e-3; // (x/a)**2 + (y/a)**2 + (z/(a/3))**2 = 1
const real radius_b = radius_a/2.0; // vol = 3/4*pi*a*a*b
const real LX = 8*radius; // width = 4*db m
const real LZ = LX*(NZ/NX);

const real gravity = -9.8; // [m/s2]
const real c_surf = 1.0; // [kg/m3] -> non-dimensionalized value

const real z_bubble_top_target = LZ - 2.0 * radius;  // target z of bubble top [m]
const real zcent = z_bubble_top_target - radius_b;  // z of initial bubbles [m]

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
  boil::plot = new PlotTECMPI();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), NX, Periodic::no());
  Grid1D gz( Range<real>(0.0, LZ), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);
  d.save("grid2048.grd");

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter gas(d), water(d);
  gas.rho   (1.977);
  gas.mu    (1.48e-5);
  gas.lambda(1.88e-9);
  water.rho   (997.);
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
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );

  Matter mixed(water, gas, & c);
  mixed.sigma(0.070);

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 10.0 * pow(gas.rho()->value()*pow(dxmin,3.0)
                        /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 10000000;
  Times time(ndt, dt);
  const real tint = 1.0e-3;  // time interval for Tecplot output
  const int nint= 1000;      //set the backup time interval
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
  enthFD.set_conv_divu_subtract(false);  // when phase change takes place

  /* phase change */
  const real latent = 1.0;
  PhaseChange pc(mdot, tpr, q, c, g, f, c, uvw, time, & mixed, latent, c_surf);
  //pc.set_rhol(1.0e+12);
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
      real aLHS = (pow(c.xc(i)/radius_a,2.0)
                  +pow(c.yc(j)/radius_a,2.0)
                  +pow((c.zc(k)-zcent)/radius_b,2.0));
      if (aLHS < 0.8) {
        c[i][j][k]=0.0;
      } else if(aLHS < 1.2) {
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
              real aLHS = pow(xxc/radius_a,2.0)
                         +pow(yyc/radius_a,2.0)
                         +pow((zzc-zcent)/radius_b,2.0);
              if (aLHS>1.0){
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

    boil::plot->plot("uvw-c-press-tpr-mdot",0,&time,&uvw,&c,&press,&tpr,&mdot);
  }

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
    w_inlet += acc * time.dt();  // inlet velocity [m/s]
    w_camera = -w_inlet;           // camera velocity [m/s]
    z_camera += w_camera * time.dt();  // camera position [m]
    boil::oout<<"main:w_inlet:time "<<time.current_time()<<" w_inlet "<<w_inlet
              <<" accl_camera "<<acc_camera<<" z_camera "<<z_camera
              <<" z_bubble_top_earth "<<z_camera+z_bubble_top_current<<"\n";
    for_m(m) {
      uvw.bc(m).modify( BndCnd( Dir::imin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::imax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::jmin(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::jmax(), BndType::wall(), 0.0, 0.0, w_inlet ) );
      uvw.bc(m).modify( BndCnd( Dir::kmax(), BndType::inlet(), 0.0, 0.0, w_inlet ) );
    }

    /* phase change */
    pc.update();
    //ns.vol_phase_change(&f);// f-> the gradient of pressure?
    mdot = 0.0;  f = 0.0;  g = 0.0;  q = 0.0;  // stop phase change
 
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

    /* advance */
    enthFD.clrold_store();  // store clrold before update; effective for restart
    conc.new_time_step();
    conc.advance();
    for_vk(c,k){
      if(c.zc(k)<0.05*LZ){
        for_vij(c,i,j){
          c[i][j][k]= 1.0;
        }
      }
    }
    /* update clr in cipcsl2 after seed, cutneck and outlet-region */
    c.bnd_update();
    c.exchange_all();
    conc.totalvol();
    conc.front_minmax();

    // area of interface
    real total_area = 0.0;
    for_vijk(c,i,j,k){
      total_area += conc.adens[i][j][k] * c.dV(i,j,k); // [m2/m3]*[m3]
    }
    boil::cart.sum_real(&total_area);

    real total_flux_CL = 0.0;
    real z_examine = LZ - 6.0*radius;
    for_vijk(tpr,i,j,k) {
      if (tpr.zn(k)<=z_examine && z_examine<tpr.zn(k+1)) {
        total_flux_CL += tpr[i][j][k]*tpr.dSz(i,j,k)*uvw[Comp::w()][i][j][k];
        // [kg/m3]*[m2]*[m/s]=[kg/s]
      }
    }
    boil::cart.sum_real(&total_flux_CL);

    real total_flux_CL1 = 0.0;
    real z_examine1 = 0.5*LZ;
    for_vijk(tpr,i,j,k) {
      if (tpr.zn(k)<=z_examine1 && z_examine1<tpr.zn(k+1)) {
        total_flux_CL1 += tpr[i][j][k]*tpr.dSz(i,j,k)*uvw[Comp::w()][i][j][k];
        // [kg/m3]*[m2]*[m/s]=[kg/s]
      }
    }
    boil::cart.sum_real(&total_flux_CL1);

    boil::oout<<"mass_transfer_rate:time= "<<time.current_time()
        <<" total_area[m2] "<<total_area
        <<" smdot_neg[kg/s] "<<pc.get_smdot_neg()
        <<" total_flux_CL[kg/s] "<<total_flux_CL
        <<" Dl[m2/s] "<<water.lambda()->value()
        <<" Sherwood1 "<<-pc.get_smdot_neg()/total_area*(2.0*radius)/water.lambda()->value()
        <<" Sherwood2 "<<total_flux_CL/total_area*(2.0*radius)/water.lambda()->value()
        <<" Sherwood3 "<<total_flux_CL1/total_area*(2.0*radius)/water.lambda()->value()
        <<"\n";

    /* solve energy equation */
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
    real cfl_current = ns.cfl_max();
    boil::oout<<"cfl= "<<time.current_time()<<" "<<cfl_current<<"\n";
    time.control_dt(cfl_current, cfl_limit, dt);

    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);

      // subtract w_inlet just for visualization, w is on the earth fixed coordinate system
      m = Comp::w();
      for_avmijk(uvw,m,i,j,k) { uvw[m][i][j][k] += w_camera; }

      boil::plot->plot("uvw-c-press-tpr-mdot",iint,&time,&uvw,&c,&press,&tpr,&mdot);

      // return uvw to the velocity reference to the moving coordinate system
      for_avmijk(uvw,m,i,j,k) { uvw[m][i][j][k] -= w_camera; }

      iint = int(time.current_time()/tint) + 1;
    }

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

    /*---------+
    |  exit ?  |
    +---------*/
    if(conc.get_zmaxft()>=LZ-0.1*radius){
       boil::oout<<"Bubble reaches to the top boundary. Exit.";
       break;
    }
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}

