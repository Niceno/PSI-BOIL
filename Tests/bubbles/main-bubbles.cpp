/*----------------------+
|                       |
|  set lcnd to 0.5 ...  |
|                       |
+----------------------*/
#include "Include/psi-boil.h"

const int DB= 128; int NB;
const int NX=  96;
const int NZ=NX*2;

const real RB = 0.002;    /* d = 4 [mm] */        
const real LX = 0.03;     /* L = 3 [cm] */
const real LZ = LX*NZ/NX; 

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
  Grid1D gx( Range<real>(-LX/2.0,     LX/2.0), NX, Periodic::yes());
  Grid1D gz( Range<real>(-LZ/6.0, 5.0*LZ/6.0), NZ, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  Times time(1000, 0.0001); /* ndt, dt */
	
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
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

  real xc[DB];
  real yc[DB];
  real zc[DB];

  boil::random_seed();
  for(int b=0; b<DB; b++) {
    
    real mind=LZ;
    int l;
    for(l=0; l<DB; l++) {

      xc[b]=boil::random_number(Range<real> (-LX/2.0+1.2*RB,    LX/2-1.2*RB));
      yc[b]=boil::random_number(Range<real> (-LX/2.0+1.2*RB,    LX/2-1.2*RB));
      zc[b]=boil::random_number(Range<real> (-LZ/6.0+1.2*RB,5.0*LZ/6-1.2*RB));

      mind=LZ;
      for(int o=0; o<b; o++) {
        real d=sqrt( (xc[b]-xc[o])*(xc[b]-xc[o]) +
                     (yc[b]-yc[o])*(yc[b]-yc[o]) +
                     (zc[b]-zc[o])*(zc[b]-zc[o]) );
        mind=boil::minr(d, mind);
      }

      if( mind > RB*2.2 ) break;
    } 
 
    NB = b+1;
    if(l == DB) {
      boil::oout << "# impossible to place more than " 
                 << b << " bubbles" << boil::endl;
      NB = b;
      break;
    }
  }
  OMS(defined bubble coordinates);

  /*----------+
  |  spheres  |
  +----------*/
  c = 0.0; /* liquid */
  for(int b=0; b<NB; b++) 
    for_vijk(c,i,j,k) {
      real x=c.xc(i)-xc[b]; real y=c.yc(j)-yc[b];  real z=c.zc(k)-zc[b];
      real dist = sqrt( x*x + y*y + z*z ); 
      if(dist < RB) c[i][j][k] = 1.0; /* gas   */
    }
  c.exchange();
  OMS(defined c);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (   0.00001);
  air  .rho   (   1.     );
  air  .sigma (   0.07   );
  water.mu    (   0.001  );
  water.rho   (1000.0    );
  water.sigma (   0.07   );

  Matter mixed(air, water, c);

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
    OPR(c.integral());
    conc.new_time_step();
    OPR(c.integral());
    conc.convection();
    OPR(c.integral());
    conc.advance();
    OPR(c.integral());
    conc.sharpen();
    OPR(c.integral());
    conc.tension(&xyz, mixed);
    OPR(c.integral());
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

    if(time.current_step() % 50 == 0) {
      boil::plot->plot(uvw,  c, "uvw,c",  time.current_step());
      boil::plot->plot(p,       "p",      time.current_step());
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
