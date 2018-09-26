#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =     2.0;
const real LY =     0.08;
const real LZ =     2.0; 

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
  Scalar c(d);
  Scalar f(d), kappa(d);

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );

  f = c.shape();


  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  Krylov * solver = new CG(d, Prec::di());
  int ndt=6280;
  int nint=628;
  real pi=acos(-1.0);
  Times time(ndt, pi/314.0); // ndt, dt

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

  for_avijk(c,i,j,k) {
    c[i][j][k] = 0.0;
    real dist = pow( (c.xc(i) - 0.5), 2 ) 
              + pow( (c.zc(k) - 0.0), 2 ); 
    dist = sqrt(dist);
    if(dist < 0.3) c[i][j][k] = 1.0;
    //if(c.xc(i)<0.6 && fabs(c.zc(k))<0.06) c[i][j][k] = 0.0;
    if(c.xc(i)<0.7 && fabs(c.zc(k))<0.06) c[i][j][k] = 0.0;
  }

  /*------------------+
  |  define a solver  |
  +------------------*/
  CIPCSL2 T(c, f, kappa, uvw, time, solver);
  T.set_nredist(1);
  T.set_itsharpen(4);
  T.totalvol();
  boil::plot->plot(c,"c", 0);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    T.new_time_step();
    //T.convection();
    T.advance();
    //T.sharpen();
    T.totalvol();

    if(boil::plot && time.current_step()%nint == 0) {
      boil::plot->plot(c,"c", time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-zalesak.cpp,v 1.4 2018/09/26 10:09:00 sato Exp $'/
+-----------------------------------------------------------------------------*/
