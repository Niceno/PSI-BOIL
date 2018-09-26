/*----------------------+
|                       |
|  set lcnd to 0.5 ...  |
|                       |
+----------------------*/
#include "Include/psi-boil.h"

/* boundary conditions */
const int NX=128;
const int NY=  4;
const int NZ= NX;

const real RB =   0.0025;           
const real LX =   RB*8.0;           
const real LZ =   LX*NZ/NX; 
const real LY =   LX*NY/NX;


/******************************************************************************/
main(int argc, char * argv[]) {

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

  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), 
             Range<real>( LY/NY,  LY/NY ),
              NY, Periodic::yes());

  Grid1D gz( Range<real>(-LZ/4.0, 3.0*LZ/4.0), 
             Range<real>( LZ/NZ,  LZ/NZ ),
              NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  Times time(200, 0.00016); /* ndt, dt */
	
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
        c[i][j][k] = 1.0; /* gas   */
      else
        c[i][j][k] = 0.0; /* liquid*/
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
  Matter air(d), water(d);
  air  .mu    (   0.00001);
  air  .rho   (   1.     );
  water.mu    (   0.001  );
  water.rho   (1000.0    );

  Matter mixed(air, water, c);
  mixed.sigma (   0.06   );

  ColorFunction conc  (c,   g, 1.0, 1.0, uvw, time, solver); 
  conc.convection_set(ConvScheme::superbee());
  conc.sharpen();
  boil::plot->plot(c, "c1", 0);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
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

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-2));
    p.exchange();
    ns.project(p);
    pr.update_rhs();

    if(time.current_step() % 5 == 0) {
      boil::plot->plot(uvw,  p, "uvw,p",  time.current_step());
      boil::plot->plot(xyz,  c, "xyz,c",  time.current_step());
    }
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw,c, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-bubble-2d.cpp,v 1.12 2012/09/13 08:13:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
