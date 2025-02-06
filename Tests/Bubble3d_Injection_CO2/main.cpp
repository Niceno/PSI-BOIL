/*----------------------------------------------------------------------------+
|  CO2 injected into the moving water                                         |
|  Experiment by Alvaro                                                       |
|  WMS: Wire Mesh Sensor
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"
//#define USE_FLOOD       // FloodFill algorithm for bubble ID, may cause crush
//#define USE_TECMPI        // Output szplt file, which only Tecplot can read.
#define OPTWMS            // Output WMS data

/* domain and cells */
const int NX= 64;       // number of cells in the X and Y directions
const real LX = 0.050;    // width of domain  = 0.05 [m]
const int NZ= NX*2;       // number of cells in the Z direction
const real LZ = LX*NZ/NX; // height of domain = 0.30 [m]

const real gravity = -9.8;// [m/s2]
const real c_surf = 1.0;  // concentration of CO2 at liquid-gas interface [-]

const real period_WMS = 1/3200.0; // period of WMS = 1.0/frame_rate [s]
//const int  nWMS = 5;     // number of WMSs
//const real zWMS[nWMS] = {0.042, 0.09, 0.14, 0.19, 0.24};//  WMSs [m]
const int  nWMS = 2;     // number of WMSs
const real zWMS[nWMS] = {0.042, 0.09};//  WMSs [m]

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
#ifdef USE_TECMPI
  boil::plot = new PlotTECMPI();
#else
  boil::plot = new PlotTEC();
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

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel [m/s]
  Scalar p  (d), f  (d), press(d,"press"); // pressure [Pa]
  Scalar c  (d,"vfl"), g  (d), kappa(d); // color function [-]
  Scalar tpr(d,"eps"), q  (d); // concentration [-]
  Scalar mdot(d,"mdot");  // phase change rate [kg/m3s]
                          // negative for dissolution of CO2 (gas) into liquid water
#ifdef USE_FLOOD
  Scalar idFlood(d,"idFlood");
#endif

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet(),0.0,0.0,"0.5*(0.3965+0.8287)-0.5*(0.8287-0.3965)*tanh(5000.0*(sqrt((x)^2+ (y)^2)-0.004))"));
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
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
  mdot = p.shape();
#ifdef USE_FLOOD
  idFlood = p.shape();
#endif

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::inlet(),"0.5*(1.0+tanh(10000.0*(sqrt((x)^2+ (y)^2)-0.004)))") );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  
  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::inlet(),"0.5*(1.0-tanh(10000.0*(sqrt((x)^2+ (y)^2)-0.004)))") );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter gas(d), water(d);
  gas.rho   (1.977);             // density of CO2 [kg/m3]
  gas.mu    (1.48e-5);           // dynamic viscosity [Pa.s]
  gas.cp    (1.0);               // must be 1.0
  gas.lambda(1.88e-9);           // D: diffusivity of salt in gas (this value doesn't make sense)
  water.rho   (997.);            // density [kg/m3]
  water.mu    (1.0e-3);          // dynamic viscosity [Pa.s]
  water.cp    (1.0);             // must be 1.0
  water.lambda(1.88e-9);         // D: diffusivity coefficient of salt in water [m2/s]
  Matter mixed(water, gas, & c); // c: volume fraction of water [-]
  mixed.sigma(0.070);            // surface tension coefficient [N/m]

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 10.0 * pow(gas.rho()->value()*pow(dxmin,3.0)
                        /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 30000;
  Times time(ndt, dt);
  const real tint = 0.001;
  const int nint= 2000; //set the backup time interval
  const real cfl_limit=0.25;

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
  conc.set_limit_color(true);

  /* enthalpy equation */
  EnthalpyFDhighPr enthFD(tpr, q, c, uvw, time, solver, & mixed ,c_surf);
  //EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, & mixed ,c_surf);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  enthFD.set_conv_divu_subtract(false);  // when phase change takes place

  /* phase change */
  const real latent = 1.0;
  PhaseChange pc(mdot, tpr, q, c, g, f, c, uvw, time, & mixed, latent, c_surf);
  pc.set_rhol(1.0e+12); // increase of liquid volume due to dissolved CO2 is negligible
  //pc.set_turbP(1.0);  // activate when turbulent eddy is introduced

  /* floodfill */
#ifdef USE_FLOOD
  Floodfill flood(c, idFlood, &uvw, time);
  flood.set_out_freq(1);
  boil::oout<<"main:Floodfill:output frequency= "<<flood.get_out_freq()<<"\n";
#endif

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
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);

    uvw.load     ("uvw",   ts);
    press.load   ("press", ts);
    c.load       ("c",     ts);
    tpr.load     ("tpr",   ts);
   
  } else {
  /*--------------------+
  |  initial condition  |
  +--------------------*/
  /* volume fraction of liquid */
  for_vijk(c,i,j,k)
    c[i][j][k] = 1.0;

  c.bnd_update();
  c.exchange_all();
  conc.init();

  /* concentration of CO2 */
  tpr.bnd_update();
  tpr.exchange_all();

#ifdef USE_FLOOD
  boil::plot->plot(uvw,c,press,tpr,mdot,idFlood, "uvw-c-press-tpr-mdot-id",0);
#else
#ifndef USE_TECMPI
  boil::plot->plot(uvw,c,press,tpr,mdot, "uvw-c-press-tpr-mdot",0);
#else
  boil::plot->plot("uvw-c-press-tpr-mdot",0,&time,&uvw,&c,&press,&tpr,&mdot);
#endif
#endif
 }

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;
  int itWMS = int(time.current_time()/period_WMS) + 1;

#if 1
    /* check mass flow rate */
    real flowrate_l = 0.0;
    real flowrate_g = 0.0;
    real area_bottom = 0.0;
    for_vijk(c,i,j,k) {
      if (approx(c.zn(k),0.0)) {
        flowrate_l += c.dSz(i,j,k) * c[i][j][k-1] * uvw[Comp::w()][i][j][k];
        flowrate_g += c.dSz(i,j,k) * (1.0-c[i][j][k-1]) * uvw[Comp::w()][i][j][k];
        area_bottom += c.dSz(i,j,k);
      }
    }
    boil::cart.sum_real(&flowrate_l);
    boil::cart.sum_real(&flowrate_g);
    boil::cart.sum_real(&area_bottom);
    boil::oout<<"main:flowrate_l= "<<flowrate_l<<" [m3/s], gas= "<<flowrate_g
              <<" [m3/s], area_bottom= "<<area_bottom<<"[m2]\n";
    boil::oout<<"main:flowrate_l= "<<1000.0*flowrate_l*3600.0<<" [l/h], gas= "
              <<1000.0*flowrate_g*60.0<<" [l/min]\n";
#endif

  /*------------+
  |  Time loop  |
  +------------*/

  for(time.start(); time.end(); time.increase()) {

    /* phase change */
    pc.update();

    /* reset body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* gravity force */ 
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }
    xyz.exchange();

    /* surface tension */
    conc.tension(&xyz, mixed);

    /* Navier-Stokes equations */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(0.001));

    p = 0.0;       // pressure difference between n+1 and n time step, dp
    multigrid.vcycle(ResRat(1e-3));
    ns.project(p); // uvw_n+1
    press += p;    // press_n+1 = press_n + dp
    press.exchange();

    /* VOF */
    enthFD.clrold_store();  // store clrold before update; effective for restart
    conc.new_time_step();
    conc.advance();
    conc.totalvol();
    conc.front_minmax();

    /* statistics */
    // area of interface
    real total_area = 0.0;
    for_vijk(c,i,j,k){
      total_area += conc.adens[i][j][k] * c.dV(i,j,k); // [m2/m3]*[m3]
    }
    boil::cart.sum_real(&total_area);
    // integrate CO2 in liquid

    real CO2_liquid = 0.0;
    for_vijk(c,i,j,k){
      CO2_liquid += c[i][j][k] * tpr[i][j][k] * c.dV(i,j,k); // [kg/m3] * [m3]
    }
    boil::cart.sum_real(&CO2_liquid);

    boil::oout<<"mass_transfer_rate:time= "<<time.current_time()<<" dt= "<<time.dt()
        <<" total_area[m2] "<<total_area<<" smdot_neg[kg/s] "<<pc.get_smdot_neg()
        <<" CO2_liquid[] "<<CO2_liquid<<"\n";

    /* solve energy equation (CO2 in this context) */
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthFD");

    /* floodfill */
#ifdef USE_FLOOD
    flood.identify_regions();
#endif

    /* dt control */
    real cfl_now = ns.cfl_max();
    boil::oout<<"cfl= "<<time.current_time()<<" "<<cfl_now<<"\n";
    time.control_dt(cfl_now, cfl_limit, dt);

    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
#ifdef USE_FLOOD
      boil::plot->plot(uvw,c,press,tpr,mdot,idFlood, "uvw-c-press-tpr-mdot-id",iint);
#else
#ifndef USE_TECMPI
      boil::plot->plot(uvw,c,press,tpr,mdot, "uvw-c-press-tpr-mdot",iint);
#else
      boil::plot->plot("uvw-c-press-tpr-mdot",iint,&time,&uvw,&c,&press,&tpr,&mdot);
#endif
#endif

#if 0
      // subtract gravity just for visualization
      m = Comp::w();
      for_vmijk(xyz,m,i,j,k) {
        xyz[m][i][j][k] -= gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
      }
      xyz.exchange();
      boil::plot->plot(xyz,c,kappa, "xyz-c-kappa",iint);
#endif
      iint = int(time.current_time()/tint) + 1;
    }

    /* output data at WMS */
#ifdef OPTWMS
    #include "outputWMS.cpp"
#endif

    if( time.current_step() % nint == 0 ) {
      uvw.save     ("uvw", time.current_step());//save store the data as the binary file '.bck'
      press.save   ("press",   time.current_step());
      c.save       ("c",   time.current_step());
      tpr.save     ("tpr",   time.current_step());
      if( boil::cart.iam()==0) {
        std::fstream output;
        output << std::setprecision(16);
        std::string name = name_file("time", ".txt", time.current_step());
        output.open(name, std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output.close();
      }
    }
    if(boil::timer.current_min() > wmin
      || time.current_step()==time.total_steps()) {
      uvw.save     ("uvw",  time.current_step());//save the data as the binary file '.bck'
      press.save   ("press",time.current_step());
      c.save       ("c",    time.current_step());
      tpr.save     ("tpr",  time.current_step());
      if( boil::cart.iam()==0) {
        std::fstream output;
        output << std::setprecision(16);
        output.open("time.txt", std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output.close();
      }
      uvw.rm  ("uvw",ts);  // rm bck files previously used for the restart
      press.rm("press",ts);
      c.rm    ("c",ts);
      tpr.rm  ("tpr",ts);
      mdot.rm ("mdot",ts);
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}

