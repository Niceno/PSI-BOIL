#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =   0.01;  // 1cm
const real LZ =   0.02;  // 2cm
const real db =   0.004; // 4mm

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
             Range<real>( LX/ 64.0, LX/ 64.0 ),
              64, Periodic::no());

  Grid1D gz( Range<real>(0,LZ), 
             Range<real>( LZ/128.0, LZ/128.0 ),
             128, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  // ok Times time(128, 0.00002); /* ndt, dt */
  Times time(1024, 0.00001); /* ndt, dt */
	
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
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
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
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  
  /*---------+
  |  circle  |
  +---------*/
  for_vijk(c,i,j,k) {
    real dist = pow( (c.xc(i)),          2 ) 
              + pow( (c.yc(j)),          2 ) 
              + pow( (c.zc(k) - LZ/2.0), 2 ); 
        dist = sqrt(dist);
        if(dist < db/2.0)
          c[i][j][k] = 1.0;
        else
          c[i][j][k] = 0.0;
   
    /* water at the bottom */
    if(c.zc(k) > 0.75 * LZ)
      c[i][j][k] = 1.0;
  }

  /*----------------------------------------+
  |  square to check surface tension force  |
  for_vijk(c,i,j,k) 
    if( fabs(c.xc(i)) < 0.15 && fabs(c.zc(k)-0.50) < 0.15 ) 
      c[i][j][k] = 0.0;
    else
      c[i][j][k] = 1.0;
  +----------------------------------------*/

  c.exchange();
  boil::plot->plot(c, "c0", 0);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (0.00001 * 1.0);
  air  .rho   (1.);
  water.mu    (0.001 * 1.0);
  water.rho   (1000.0);

  Matter mixed(air, water, c);

  ColorFunction conc  (c,   g, 1.0, 1.0, uvw, time, solver); 
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

  /*----------+
  |  gravity  |
  +----------*/
  const Comp m = Comp::w();
  for(int i=1; i<xyz.ni(m)-1; i++)
    for(int j=1; j<xyz.nj(m)-1; j++)
      for(int k=1; k<xyz.nk(m)-1; k++)
        xyz[m][i][j][k] = -9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.convection();
    ns.solve(ResRat(0.0001));

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-2));
    p.exchange();
    ns.project(p);
    pr.update_rhs();

    /* advance */
    conc.new_time_step();
    conc.convection();
    conc.advance();
    conc.sharpen();
    conc.tension(&xyz, mixed);
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] -= 9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    ns.cfl_max();

    /*------------+
    |  save/plot  |
    +------------*/
    if(time.current_step() % 16 == 0) {
      boil::plot->plot(uvw,  c, "uvw,c",  time.current_step());
      boil::plot->plot(p,    "p",    time.current_step());
      uvw.  save("uvw",  time.current_step());
      p.    save("p",    time.current_step());
      ns.   save("ns",   time.current_step());
      pr.   save("pr",   time.current_step());
      conc. save("conc", time.current_step());
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
/*-----------------------------------------------------------------------------+
 '$Id: main-drop.cpp,v 1.18 2012/09/13 08:13:29 niceno Exp $'/
+-----------------------------------------------------------------------------*/
