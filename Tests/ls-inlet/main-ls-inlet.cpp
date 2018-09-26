#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =    10.0;
const real LY =     0.5;
const real LZ =     0.5; 
const real mu =   0.1;
const real dt =   0.005; 

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
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), 200, Periodic:: no() );
  Grid1D gy( Range<real>(-0.5*LY,0.5*LY),  10, Periodic:: no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gy);
	
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

  Times time(5000, 0.001); // ndt, dt

  std::cout << phi.bc() << std::endl;
  OMS(adding b.c.);
  phi.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 1.0 ) );
  phi.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );

  /*------------------+
  |  define a solver  |
  +------------------*/
  ColorFunction T(phi, f, 1.0, 1.0, uvw, time, NULL); 

  /*------------------------+
  |  divergence-free field  |
  +------------------------*/
  const Comp m = Comp::u();
  for(int i=0; i<uvw.ni(); i++)
    for(int j=0; j<uvw.nj(); j++)
      for(int k=0; k<uvw.nk(); k++) 
        uvw[m][i][j][k] =  1.0;          
  
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

    if(boil::plot && time.current_step()%500 == 0) {
      boil::plot->plot(uvw,"uvw", time.current_step());
      boil::plot->plot(phi,"phi", time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(phi, "test", 0);
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-ls-inlet.cpp,v 1.10 2012/09/13 08:13:33 niceno Exp $'/
+-----------------------------------------------------------------------------*/
