#include "Include/psi-boil.h"

const int NX =  64;
const int NZ = NX * 4;
const Range<real> DB(0.0005, 0.001);    /* diameter = 0.5 - 1 mm */
const real        DB_AVG = (DB.first()+DB.last()) / 2.0;
const real        LX = DB_AVG *  8.0;
const real        LZ = LX * NZ / NX; 

/* time steps */
const int ndt  = 10240; /* number of time steps */
const int sint =    64;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0,     LX/2.0), NX, Periodic::no());
  Grid1D gz( Range<real>(-LZ/8.0, 7.0*LZ/8.0), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air.rho  (1.205);
  air.mu   (1.82e-5);   
  water.rho(998.2);
  water.mu (1.0e-3); 

  real responce_time = 1.205 * DB_AVG * DB_AVG / (36.0*1.0e-3);
  real dt = responce_time / 24.0;     
  OPR(dt);
  Times time(ndt, 0.00005); /* ndt, dt */

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d), c_m(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d), step(d), sflag(d); // concentration
  Scalar press(d);
	
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

  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  for_m(m) {
    c_m.bc(m).add( BndCnd( Dir::imin(), BndType::neumann() ) );
    c_m.bc(m).add( BndCnd( Dir::imax(), BndType::neumann() ) );
    c_m.bc(m).add( BndCnd( Dir::jmin(), BndType::neumann() ) );
    c_m.bc(m).add( BndCnd( Dir::jmax(), BndType::neumann() ) );
    c_m.bc(m).add( BndCnd( Dir::kmin(), BndType::neumann() ) );
    c_m.bc(m).add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  }

  Matter mixed(water, air, c);
  mixed.sigma(0.072);
  Dispersed disp(c, /* g, */ uvw, time, &mixed); 

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns(uvw, xyz, time, solver, &mixed);  
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());

  Pressure pr(p, f, uvw, time, solver, &mixed); 

  AC multigrid( &pr );

  const real D0 = 0.0005;
  const real D1 = 0.00075;
  const real D2 = 0.001; 
  const real D3 = 0.00125; 
  disp.add( Particle( Position(  0.0,  0.0, 0.0), Diameter(D3) ) );
  const real offst = D3 * 0.3535; /* sqrt(2)/4 */
  disp.add( Particle( Position( -offst, -offst,     D3), Diameter(D0) ) );
  disp.add( Particle( Position(  offst, -offst, 2.0*D3), Diameter(D0) ) );
  disp.add( Particle( Position(  offst,  offst, 3.0*D3), Diameter(D0) ) );
  disp.add( Particle( Position( -offst,  offst, 4.0*D3), Diameter(D0) ) );
    
  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /*------------------------------+
    |  add new bubble from sparger  |
    if( (time.current_step()-1) % bint == 0) {
      disp.add( Particle (Position(boil::random_number(SP),
                                   boil::random_number(SP),
                                   0.0), 
                          Diameter(boil::random_number(DB))) );
    }
    +------------------------------*/

    for_m(m) 
      for_vmijk(xyz,m,i,j,k) 
        xyz[m][i][j][k] = 0.0;

    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = -9.81 * xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }

    disp.advance(& xyz); /* will add forces on bubbles */

    if( ns.cfl_max() > 1.0 ) {
      boil::plot->plot(uvw,c, press, "uvw-c-press",time.current_step());
      boil::oout << "cfl too high, stopping!" << boil::endl;
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
 
    ns.grad(press);
    ns.convection();
    ns.solve(ResRat(1e-2));   

    p = 0.0;
    p.exchange();
    multigrid.vcycle(ResRat(1e-2));  
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();
    pr.update_rhs();

    /* time step control */
    time.control_dt(ns.cfl_max(), 0.33, 1.0);
    OPR(dt);
    OPR(time.dt());

    if(time.current_step() % sint == 0) {
      boil::plot->plot(uvw,c, press, "uvw-c-press",time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();
}	
