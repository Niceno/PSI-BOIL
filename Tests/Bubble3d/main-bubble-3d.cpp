/*----------------------+
|                       |
|  set lcnd to 0.5 ...  |
|                       |
+----------------------*/
#include "Include/psi-boil.h"

/* boundary conditions */
const int NX= 96;
const int NZ= NX*2;

const real RB = 0.0025;           
const real LX = RB*6.0;           
const real LZ = LX*NZ/NX; 

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC(AsNodes::yes());

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 
             Range<real>( LX/NX,  LX/NX ),
              NX, Periodic::no());

  Grid1D gz( Range<real>(-LZ/6.0, 5.0*LZ/6.0), 
             Range<real>( LZ/NZ,  LZ/NZ ),
              NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  Times time(10, 0.00004); /* ndt, dt */
	
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

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet(), 0.0, 0.0, 0.5 ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  
  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 0.0 ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /*---------+
  |  circle  |
  +---------*/
  for_vijk(c,i,j,k) {
    real dist = sqrt( c.xc(i)*c.xc(i) + c.yc(j)*c.yc(j) + c.zc(k)*c.zc(k) ); 
      if(dist < RB)
        c[i][j][k] = 1.0; /* gas   */
      else
        c[i][j][k] = 0.0; /* liquid*/
  }

  c.exchange();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (   0.00001);
  air  .rho   (   1.2    );
  water.mu    (   0.001  );
  water.rho   ( 999.0    );

  Matter mixed(air, water, c);
  mixed.sigma (   0.07   );

  ColorFunction conc  (c,   g, 1.0, 1.0, uvw, time, solver); 
  conc.convection_set(ConvScheme::superbee());
  conc.sharpen();

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());

  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  const Comp m = Comp::w();
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k] = 0.5;
  boil::plot->plot(uvw, c, "initial", 0);

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
   
    real tot=0.0;
    for_vijk(c,i,j,k)
      tot += c.dV(i,j,k) * c[i][j][k];
    OPR(tot);
 
    const Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] -= 9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    ns.cfl_max();

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();

    ns.convection();
    ns.solve(ResRat(0.01));
    // automatic: ns.get_src(&f); ???

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-2));
    p.exchange();
    ns.project(p);
    // automatic: ns.get_src(&f); ???

    if(time.current_step() % 10 == 0) {
      boil::plot->plot(uvw,  c, "uvw,c",  time.current_step()/10);
      boil::plot->plot(p,       "p",      time.current_step()/10);
    }

    /*---------+
    |  exit ?  |
    +---------*/
    std::ifstream infile;
    infile.open ("stop.now", std::ifstream::in);
    if(infile.good()) exit(0);
    infile.close();
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw,c, "test", 0);
}	
