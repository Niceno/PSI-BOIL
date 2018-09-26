#include "Include/psi-boil.h"

/* boundary conditions */
const real LX_IN =  0.1; 
const real LX    = 32.0; 
const real LY    =  2.0;
const real LZ    =  4.0;
const real mu    =  0.0003; // 0.001, 0.0005 was OK; target: 0.0003

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx_in(Range<real>( 0.0, LX_IN ), 2, Periodic::yes());

  Grid1D gx(Range<real>( 0.0, LX ), 
            Range<real>( LX/256.0, LX/256.0 ),
           256, Periodic::no());

  Grid1D gy(Range<real>( -LY/2.0,  LY/2.0), 
            Range<real>( LY/ 64.0, LY/ 64.0 ),
            32, Periodic::no());

  Grid1D gz(Range<real>( 0.0, LZ ), 
            Range<real>( LZ/64.0, LZ/64.0 ),
            64, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d_in(gx_in, gy, gz, "inlet_dom", Decompose::no());
  Domain d   (gx,    gy, gz);
	
  Krylov * solver = new CG(d, Prec::di());

  Times time(200, 0.03); /* ndt, dt 0.015 */

  int ne = 100;
  real tt = 0.1; // turbulent time scale
  real tl = 0.1; // turbulent length scale
  RandomFlow rf( ne, tt, tl );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu(mu);
	
  Vector uvw_in(d_in);   // inlet velocity
  Vector uvw(d), xyz(d); // computed velocity
  Scalar p(d),   f(d);

  real uvw_bulk[3] = {1.0, 0.0, 0.0};

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 1.0, 0.0, 0.0 ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  Momentum ns(uvw, xyz,      time, solver, & fluid);
  Pressure pr(p,   f,   uvw, time, solver, & fluid);

  AC sol( &pr );

  boil::plot = new PlotTEC();

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    /* inlet plane */
    //uvw_in.randomize(rf, 0.2, 0.0);
    uvw_in.randomize(rf, 0.5, time.current_time());
    uvw_in += uvw_bulk;
    uvw.copy(uvw_in, Dir::imin(),1,1);

    //boil::plot->plot(uvw_in,"uvw_in", time.current_step());


    /* true domain */
    ns.cfl_max();
    ns.new_time_step();

    ns.solve(ResRat(0.0001));

    for(int i=0; i<p.ni(); i++)
      for(int j=0; j<p.nj(); j++)
        for(int k=0; k<p.nk(); k++)
          p[i][j][k] = 0.0;
    
    OMS(about to  start cycle);
    sol.vcycle(ResRat(0.01));
    p.exchange();
    ns.project(p);

    /*-------+
    |  save  |
    +-------*/
    if( time.current_step() % 200 == 0 ) {
      boil::plot->plot(uvw,"uvw", time.current_step());
      boil::plot->plot(p,  "p",   time.current_step());
      //uvw.save("uvw", time.current_step());
      //p.save  ("p",   time.current_step());
    }
  }

  boil::plot->plot(uvw,"uvw", time.current_step()-1);
  boil::plot->plot(p  ,"p",   time.current_step()-1);

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(p, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-random-plane.cpp,v 1.16 2011/05/25 11:14:02 niceno Exp $'/
+-----------------------------------------------------------------------------*/
