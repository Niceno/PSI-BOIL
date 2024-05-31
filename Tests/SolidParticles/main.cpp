/*----------------------------------------------------------------------------+
|  Rising air bubble in stagnant oil                                          |
|  Publication: Y.Sato, B.Niceno, Int J Numer Meth Fluids, 70 (2012) 441-467  |
|  Case (d)-de/20 in Fig.23                                                   |
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"
#define USE_VOF

void update_step(const Scalar & c, Scalar & step, Scalar & sflag);

const int NX= 128;
const int NZ= NX*1;
const real radius = 0.5 * 0.0261;  // bubble diameter = 0.0261 m
const real LX = radius*12.0;
const real LZ = LX*NZ/NX; 

const real gravity = -9.8;

const real p_den = 1.0e+3;  // particle density
const real p_dia = 1.0e-4;  // particle diameter

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

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity and body force
  Scalar p  (d,"dp"), f  (d), press(d,"pressure"); // pressure
  Scalar c  (d,"color"), g  (d), kappa(d); // concentration
#ifndef USE_VOF
  Scalar step(d), sflag(d);
#endif

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
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
  
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );

#ifndef USE_VOF
  step = c.shape();
  sflag = c.shape();
#endif

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), oil(d);
  air.rho   (1.205);
  air.mu    (1.82e-5);
  oil.rho   (1350.);
  oil.mu    (0.7745); //(d)
#ifdef USE_VOF
  Matter mixed(oil, air, & c);
#else
  Matter mixed(oil, air, & step);
#endif
  mixed.sigma(0.078);

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 5.0 * pow(air.rho()->value()*pow(dxmin,3.0)
                        /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 10000;
  Times time(ndt, dt);
  const real tint = 0.01;
  const int nint= 100;
  const real cfl_limit=0.25;

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  /*  linear solver */
  Krylov * solver = new CG(d, Prec::ic2());
  /* Navier-Stokes */
  Momentum ns( uvw, xyz, time, solver, &mixed);
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.max_cycles(10);
  multigrid.min_cycles(4);

  /* VOF or CIPCSL2 */
#ifdef USE_VOF
  VOF conc (c,   g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc (c,   g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
#endif

  /* SolidParticle */
  SolidParticle solidp(uvw, & time, & mixed, 0.0, 0.0, gravity, & c, & press);

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

    uvw.load     ("uvw",     ts);
    press.load   ("press",   ts);
    c.load       ("c",       ts);
    solidp.load  ("solidp",ts);

  } else {
    /*--------------------+
    |  initial condition  |
    +--------------------*/

    /* color function */
    // initialize (not mandatory)
    c = 0.0;

    // define bubble center
    const real xcent =0.0;
    const real ycent =0.0;
    const real zcent =LZ*0.1;
    // define sphere (inside = 1, outside = 0)
    boil::setup_sphere(c, radius, xcent, ycent, zcent, 8);

    // reverse
    for_vijk(c,i,j,k)
      c[i][j][k] = 1.0 - c[i][j][k];
    c.exchange_all();
    c.bnd_update();

    /* solidparticle */
    // pattern 1: define in each decomposed domain
    int np_local = 0;
    for_vijk(c,i,j,k) {
      if (c[i][j][k]< 0.5) {
         if (i%2==0 && j%2==0 && k%2==0) {
           solidp.add_local(c.xc(i), c.yc(j), c.zc(k), & p_dia, & p_den);
	   np_local++;
         }
      }
    }
    solidp.exchange();  // Never forget this!
    std::cout<<"main:np_local= "<<np_local<<"\n";
    boil::cart.sum_int(&np_local);
    boil::oout<<"main:np_global= "<<np_local<<"\n";

    // pattern 2: define by coordinates
    for (int i = 0; i < 32; i++) {
      for (int j = 0; j < 32; j++) {
        real xx = 2.0 * dxmin * i;
        real yy = 2.0 * dxmin * j;
        solidp.add_global(xx, yy, 0.5*LZ, & p_dia, & p_den);
      }
    }

    solidp.init();

    boil::plot->plot(uvw,c,press, "uvw-c-press",0);
    boil::plot->plot(solidp, "particles",0);
    conc.front_minmax();
  }
#if 1
#ifndef USE_VOF
  boil::update_step(c, step, sflag);
#endif

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  Time loop  |
  +------------*/

  for(time.start(); time.end(); time.increase()) {

    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;
 
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }
    xyz.exchange();


#ifdef USE_VOF
    conc.tension(&xyz, mixed);
#else
    conc.tension(&xyz, mixed, step);
#endif

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
    conc.new_time_step();
    conc.advance();
#ifndef USE_VOF
    boil::update_step(c, step, sflag);
#endif
    conc.totalvol();
    conc.front_minmax();

    /* update solidp */
    boil::oout<<"main:before advance\n";
    solidp.advance();
    boil::oout<<"main:after advance\n";

    /* dt control */
    time.control_dt(ns.cfl_max(), cfl_limit, dt);

    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      //boil::plot->plot(uvw,c,press, "uvw-c-press",iint,&time);
      //boil::plot->plot(solidp, "solidp",iint,&time);
      boil::plot->plot(uvw,c,press, "uvw-c-press",iint);
      boil::plot->plot(solidp, "particles",iint);

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
    if( time.current_step() % nint == 0 ) {
      // output variables (*.bck) for restart
      uvw.save     ("uvw",     time.current_step());
      press.save   ("press",   time.current_step());
      c.save       ("c",       time.current_step());
      solidp.save("solidp",time.current_step());
      // output time-*.txt file
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

    /*---------+
    |  exit ?  |
    +---------*/
    if(conc.get_zmaxft()>=LZ*0.95){
       boil::oout<<"Bubble reaches to the top boundary. Exit.";
       break;
       //exit(0);
    }
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
#endif
}

