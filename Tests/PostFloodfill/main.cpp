/*----------------------------------------------------------------------------+
|  Injection of gas into stagnant liquid                                      |
|  K. Hayashi, A. Tomiyama, Interface Tracking Simulation of Mass Transfer    |
|  from a Dissolving Bubble, The Journal of Computational Multiphase Flows,   |
|  3 (2011) 247-261. 10.1260/1757-482x.3.4.247                                |
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"

const int NX= 64;
const int NZ= NX*3;
const real RB = 0.5 *0.01;  // diameter = 0.005 m
const real LX = 0.048; // width = 0.048 m
const real LZ = LX*NZ/NX; // height = 0.24 m

const real gravity = -9.8;
const real c_surf = 1.452;

const real Rp_in  = 0.003;  // inner pipe radius
const real Rp_out = 0.004;  // inner pipe radius
const real H_pipe = 0.010;   // hight of pipe
const real wgas   = 0.51283;  // gas inlet velocity
const real wwater = 0.36618;  // water inlet velocity

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
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 
             Range<real>( LX/NX,  LX/NX ),
              NX, Periodic::no());

  Grid1D gz( Range<real>(0.0, LZ), 
             Range<real>( LZ/NZ,  LZ/NZ ),
              NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air.rho   (1.205);
  air.mu    (1.48e-5);
  air.cp    (1.0);
  air.lambda(1.33e-9);
  water.rho   (997.);
  water.mu    (1.0e-3); //dynamic viscosity [Pa.s]
  water.cp    (1.0);
  water.lambda(1.88e-9);// D diffusivity coefficient of salt in water [m2/s]

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d), press(d,"press"); // pressure
  Scalar c  (d,"vfl"), g  (d), kappa(d); // color function
  Scalar tpr(d,"eps"), q  (d); // concentration
  Scalar idFlood(d,"idFlood");

  /*-----------------------------+ 
  |  insert boundary conditions  |
 +-----------------------------*/

  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    //uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet(),0.0,0.0,0.2));
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet(),0.0,0.0,"0.5*(0.87902-0.14665*tanh(10000.0*(sqrt((x)^2+ (y)^2)-0.003)))"));
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  }

/*
const real RADIUS = 0.004;
const real VELOCITY = 0.332;
const real CENTER_X = 0.0;
const real CENTER_Y = 0.0;

for_m(m) {
  for(int i = 1; i <= NX; ++i) {
    for(int j = 1; j <= NX; ++j) {
//   for_vmij(uvw,m,i,j) 
      real x = gx.xc(i);
      real y = gx.xc(j);
      real distance = sqrt((x - CENTER_X)*(x - CENTER_X) + (y - CENTER_Y)*(y - CENTER_Y));
      if(distance <= RADIUS) {
        uvw.bc(m).add(BndCnd(Range<int>(i,i),Range<int>(j,j),Dir::kmin(), BndType::inlet(), 0.0, 0.0, VELOCITY));
      } else {
        uvw.bc(m).add(BndCnd(Range<int>(i,i),Range<int>(j,j),Dir::kmin(), BndType::wall()));
      }
  }
}
  uvw.bc(m).add(BndCnd(Dir::kmax(), BndType::outlet()));
  uvw.bc(m).add(BndCnd(Dir::jmin(), BndType::wall()));
  uvw.bc(m).add(BndCnd(Dir::jmax(), BndType::wall()));
  uvw.bc(m).add(BndCnd(Dir::imin(), BndType::wall()));
  uvw.bc(m).add(BndCnd(Dir::imax(), BndType::wall()));
}
*/
  
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
  idFlood = p.shape(); 

 /*
 for(int i = 1; i <= NX; ++i) {
   for(int j = 1; j <= NX; ++j) {
     real x = gx.xc(i);
     real y = gx.xc(j);
     real distance = sqrt((x - CENTER_X)*(x - CENTER_X) + (y - CENTER_Y)*(y - CENTER_Y));
     if(distance <= RADIUS) {
     c.bc().add(BndCnd(Range<int>(i,i),Range<int>(j,j),Dir::kmin(), BndType::dirichlet(), 0.0));
     boil::oout << "Setting boundary conditions for i=" << i << ", j=" << j << ", distance=" << distance << boil::endl;
     } else {
       c.bc().add(BndCnd(Range<int>(i,i),Range<int>(j,j),Dir::kmin(), BndType::wall()));
     }
  }
}
 boil::oout<<"vlf is defined successfully\n";

 c.bc().add(BndCnd(Dir::kmax(), BndType::outlet()));
 c.bc().add(BndCnd(Dir::jmin(), BndType::wall()));
 c.bc().add(BndCnd(Dir::jmax(), BndType::wall()));
 c.bc().add(BndCnd(Dir::imin(), BndType::wall()));
 c.bc().add(BndCnd(Dir::imax(), BndType::wall()));
*/
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  //c.bc().add( BndCnd( Dir::kmin(), BndType::inlet(),1.0) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::inlet(),"0.5*(1.0+tanh(10000.0*(sqrt((x)^2+ (y)^2)-0.003)))") );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  
  tpr.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(), 0.0 ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(), 0.0 ) );
  //tpr.bc().add( BndCnd( Dir::kmin(), BndType::inlet(),0.0) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::inlet(),"0.5*(1.0-tanh(10000.0*(sqrt((x)^2+ (y)^2)-0.003)))") );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );

  Matter mixed(water, air, & c);
  mixed.sigma(0.078);

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 5.0 * pow(air.rho()->value()*pow(dxmin,3.0)
                        /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 50000;
  Times time(ndt, dt);
  const real tint = 0.001;
  const int nint= 5000; //set the backup time interval
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

  /* Pathline */
//  Pathline pathline(uvw, & time, & c, & tpr);

  /* enthalpy equation */
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, & mixed ,c_surf);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  /* floodfill */
  Floodfill flood(c, idFlood, &uvw, time);
  flood.set_out_freq(1);
  boil::oout<<"main:Floodfill:output frequency= "<<flood.get_out_freq()<<"\n";

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
    //pathline.load("pathline",ts);
    tpr.load     ("tpr",   ts);
   
  } else {
  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vijk(c,i,j,k)
    c[i][j][k] = 1.0;

  c.bnd_update();
  c.exchange_all();

  conc.init();
  conc.front_minmax();
  
  tpr.bnd_update();
  tpr.exchange_all();

  boil::plot->plot(uvw,c,press,tpr,idFlood, "uvw-c-press-tpr-id",0);
 }

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  Time loop  |
  +------------*/

  for(time.start(); time.end(); time.increase()) {

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
    exit;
#endif

    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;
 
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }
    xyz.exchange();

    /* surface tension */
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

#if 0
    // crude approach!
    for_m(m){
      for_vmijk(uvw,m,i,j,k) {
        real rtmp = sqrt(pow(uvw.xc(m,i),2)+pow(uvw.yc(m,j),2));
        if (rtmp<Rp_in && uvw.zc(m,k)<H_pipe) {
          if (m==Comp::u()||m==Comp::v()) uvw[m][i][j][k]=0.0;
          if (m==Comp::w()) uvw[m][i][j][k]=wgas;
        }
      }
    }
#endif
    uvw.exchange_all();

    /* advance */
    conc.new_time_step();
    conc.advance();
    conc.totalvol();
    conc.front_minmax();

    /* delete bubble */
    real z_delete = LZ - 3.0 * RB;
    real cmin_delete = 1.0;
    for_vk(c,k) {
      if (c.zn(k)<=z_delete && z_delete<c.zn(k+1)) {
        for_vij(c,i,j) {
          if (cmin_delete > c[i][j][k]) cmin_delete = c[i][j][k];
        }
      }
    }
    boil::cart.min_real(& cmin_delete);
    boil::oout<<"main:cmin_delete "<<time.current_time()<<" "<<cmin_delete<<"\n";
    if (cmin_delete > 0.95) {
      for_vijk(c,i,j,k) {
        if(c.zc(k)>z_delete) c[i][j][k]=1.0;
      }
    }
    c.bnd_update();
    c.exchange_all();

    /* solve energy equation */
    enthFD.discretize();     // mu_t
    enthFD.new_time_step();  // mu_t
    enthFD.solve(ResRat(1e-16),"enthFD");

    /* floodfill */
    flood.identify_regions();

    /* dt control */
    time.control_dt(ns.cfl_max(), cfl_limit, dt);

    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,c,press,tpr,idFlood, "uvw-c-press-tpr-id",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    if( time.current_step() % nint == 0 ) {
      uvw.save     ("uvw", time.current_step());//save store the data as the binary file '.bck'
      press.save   ("press",   time.current_step());
      c.save       ("c",   time.current_step());
      tpr.save     ("eps",   time.current_step());
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
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }

    /*---------+
    |  exit ?  |
    +---------*/
#if 0 
   if(conc.get_zmaxft()>=LZ*0.99){
       boil::oout<<"Bubble reaches to the top boundary. Exit.";
       break;
    }
#endif
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}

