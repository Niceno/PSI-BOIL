#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =     2.0;
const real LY =     0.1;
const real LZ =     2.0; 
const real mu =   0.1;
const real dt =   0.005; 

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*------------------+
  |  plotting format  |
  +------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), 100, Periodic::no() );
  Grid1D gy( Range<real>(-0.5*LY,0.5*LY),   4, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);
	
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

  phi.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  phi.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  Krylov * solver = new CG(d, Prec::di());
  Times time(3140, 0.002); // ndt, dt

  /*------------------+
  |  define a solver  |
  +------------------*/
  ColorCIP T(phi, f, 1.0, 1.0, uvw, time, solver);

  /*------------------------+
  |  divergence-free field  |
  +------------------------*/
  for(int i=0; i<uvw.ni(); i++)
    for(int j=0; j<uvw.nj(); j++)
      for(int k=0; k<uvw.nk(); k++) {
        uvw[Comp::u()][i][j][k] =  uvw.zc(Comp::u(),k);
        uvw[Comp::w()][i][j][k] = -uvw.xc(Comp::w(),i);
      }
  uvw.exchange(); //set periodic boundary condition.

  for_avijk(phi,i,j,k) {
    phi[i][j][k] = 0.0;
    real dist = pow( (phi.xc(i) - 0.5), 2 ) 
              + pow( (phi.zc(k) - 0.0), 2 ); 
    dist = sqrt(dist);
    if(dist < 0.3) phi[i][j][k] = 1.0;
    if(phi.xc(i)<0.6 && fabs(phi.zc(k))<0.06) phi[i][j][k] = 0.0;
  }
      
  //T.sharpen();
  T.totalvol();
  boil::plot->plot(phi,"phi", 0);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    T.new_time_step(); /* preparations for new time step */
    T.convection();    /* computes convective term for current time step;
                          important only for implicit convection treatment */
    T.advance();       /* integrates in time explicit schemes */
    T.totalvol();

    if(boil::plot && time.current_step()%157 == 0) {
    //if(boil::plot && time.current_step()%1570 == 0) {
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
 '$Id: main-cip.cpp,v 1.1 2009/07/31 12:38:56 sato Exp $'/
+-----------------------------------------------------------------------------*/
