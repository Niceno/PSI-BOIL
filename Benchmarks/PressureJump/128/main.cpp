#include "Include/psi-boil.h"

/* boundary conditions */
const int NX=128;
const int NY=  4;

const real LX =   0.0075;  // 7.5 mm
const real LY =   LX*NY/NX;
const real RB =   0.002;   // 2mm

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
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 
             Range<real>( LX/NX,  LX/NX ),
              NX, Periodic::no());

  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), 
             Range<real>( LY/NY,  LY/NY ),
              NY, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);

  Times time(32, 0.0002); /* ndt, dt */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d); // concentration
  Scalar press(d);

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
  
  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /*---------+
  |  circle  |
  +---------*/
  for_vijk(c,i,j,k) {
    real dist = sqrt( c.xc(i)*c.xc(i) + c.zc(k)*c.zc(k) ); 
      if(dist < RB)
        c[i][j][k] = 0.0; /* liqud */
      else
        c[i][j][k] = 1.0; /* gas   */
  }

  /*----------------------------------------+
  |  square to check surface tension force  |
  for_vijk(c,i,j,k) 
    if( fabs(c.xc(i)) < RB && fabs(c.zc(k)-LZ/2.0) < RB ) 
      c[i][j][k] = 1.0;
    else
      c[i][j][k] = 0.0;
  +----------------------------------------*/

  c.exchange();
  boil::plot->plot(c, "c0", 0);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter gas(d), liquid(d);
  gas   .mu    (  0.00001);
  gas   .rho   (  1.1768 );
  gas   .sigma (  0.00   );
  liquid.mu    (  0.0012 );
  liquid.rho   (787.88   );
  liquid.sigma (  0.02361);

  Matter mixed(gas, liquid, c);

  LevelSet      conc  (c,   g, 1.0, 1.0, uvw, time, solver); 
  conc.convection_set(ConvScheme::superbee());
  conc.sharpen();
  boil::plot->plot(c, "c1", 0);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());

  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    /* advance */
    conc.new_time_step();
    conc.convection();
    conc.advance();
    conc.sharpen();
    conc.tension(&xyz, mixed);
//  boil::plot->plot(xyz,  c, "xyz,c",  time.current_step());

    ns.cfl_max();

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();

    // NEW: SIMPLE ITERATIONS
    for(int in=0; in<6; in++) {
      ns.grad(press);
      ns.convection();
      ns.solve(ResRat(0.001));

      for_vijk(p,i,j,k) p[i][j][k] = 0.0;
      p.exchange();

      multigrid.vcycle(ResRat(1e-2));
      p.exchange();
      ns.project(p);

      for_vijk(p,i,j,k)
        press[i][j][k] += p[i][j][k];
      press.exchange();
      if( pr.update_rhs() < 1.e-12) break;
    }

    if(time.current_step() % 8 == 0) {
      boil::plot->plot(uvw,  p, "uvw,p",  time.current_step());
      boil::plot->plot(xyz,  c, "xyz,c",  time.current_step());
      boil::plot->plot(press,"press",time.current_step());
    }

    /*--------------------+
    |  compute pressures  |
    +--------------------*/
    real p0 = 0.0;
    real p1 = 0.0;
    real v0 = 0;
    real v1 = 0;
    for_vijk(press,i,j,k) {
      if(press[i][j][k] < 0.5) {
        p0 += press[i][j][k] * press.dV(i,j,k);
        v0 += press.dV(i,j,k);
      } else {
        p1 += press[i][j][k] * press.dV(i,j,k);
        v1 += press.dV(i,j,k);
      }
    }
    p0 /= v0;
    p1 /= v1;
    OPR(v0+v1);
    OPR(p0-p1);

  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw,c, "test", 0);
}	

/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.10 2011/05/25 11:09:47 niceno Exp $'/
+-----------------------------------------------------------------------------*/
