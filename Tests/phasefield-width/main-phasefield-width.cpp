#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =   1.0;
const real LY =   1.0;
const real LZ =   1.0; 
const int  NX = 128;
const int  NY = 128;
const int  NZ = 128;
const real mu =   0.1;
const real dt =   0.005; 

#define DIR 0 /* 0, 1 or 2 */

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*------------------+
  |  plotting format  |
  +------------------*/
  boil::plot = new PlotGMV();

  /*----------+
  |  grid(s)  |
  +----------*/
  #if DIR==0
    Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no() );
    Grid1D gy( Range<real>(-0.5*LY,0.5*LY), NY, Periodic::yes() );
    Grid1D gz( Range<real>(-0.5*LZ,0.5*LZ), NZ, Periodic::yes() );
  #elif DIR==1
    Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::yes() );
    Grid1D gy( Range<real>(-0.5*LY,0.5*LY), NY, Periodic::no() );
    Grid1D gz( Range<real>(-0.5*LZ,0.5*LZ), NZ, Periodic::yes() );
  #else
    Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::yes() );
    Grid1D gy( Range<real>(-0.5*LY,0.5*LY), NY, Periodic::yes() );
    Grid1D gz( Range<real>(-0.5*LZ,0.5*LZ), NZ, Periodic::no() );
  #endif

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // velocity
  Scalar phi(d); // level set
  Scalar f(d);   // ???

  #if DIR==0
    phi.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
    phi.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
    phi.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  #elif DIR==1
    phi.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
    phi.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
    phi.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  #else
    phi.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    phi.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
    phi.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  #endif

  Krylov * solver = new CG(d, Prec::di());
  Times time(1, 0.002); // ndt, dt

  /*------------------+
  |  define a solver  |
  +------------------*/
  PhaseField T(phi, f, 1.0, 1.0, uvw, time, solver); // 0.0025 was OK
  T.convection_set(ConvScheme::superbee());

  /*----------------+
  |  initial field  |
  +----------------*/
  for_avijk(phi,i,j,k) {
    phi[i][j][k] = 0.0;
    #if DIR==0
      if(phi.xc(i) < 0.0) phi[i][j][k] = 1.0;
    #elif DIR==1
      if(phi.yc(j) < 0.0) phi[i][j][k] = 1.0;
    #else
      if(phi.zc(k) < 0.0) phi[i][j][k] = 1.0;
    #endif
  }
      
  T.sharpen();
  boil::plot->plot(phi,"phi", 0);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    T.new_time_step();
    T.convection();
    T.advance();
    T.sharpen();

    if(boil::plot && time.current_step()%20 == 0) {
      boil::plot->plot(uvw,"uvw", time.current_step());
      boil::plot->plot(phi,"phi", time.current_step());
    }
  }
  //T.local();

  #if DIR==0
    Rack ri("phi-x", d, Range<int>(1,NX), NY/2, NZ/2);
    ri.print(phi);
  #elif DIR==1
    Rack rj("phi-y", d, NX/2, Range<int>(1,NY), NZ/2);
    rj.print(phi);
  #else
    Rack rk("phi-z", d, NX/2, NY/2, Range<int>(1,NZ));
    rk.print(phi);
  #endif

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw,"uvw", time.current_step()-1);
  boil::plot->plot(phi,"phi", time.current_step()-1);
}	

/*-----------------------------------------------------------------------------+
 '$Id: main-phasefield-width.cpp,v 1.1 2009/12/29 10:26:38 l_niceno Exp $'/
+-----------------------------------------------------------------------------*/
