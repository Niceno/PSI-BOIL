#include "Include/psi-boil.h"
#define USE_VOF

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

  c.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(),0.0 ) );

  f = c.shape();


  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::inlet() ) );
  }

  Krylov * solver = new CG(d, Prec::di());
  //int ndt=6280;
  //int nint=628;
  int ndt=6280*4;
  int nint=628;
  //int ndt=680;
  //int nint=20;
  real pi=acos(-1.0);
  //real dt=pi/314.0; // original
  real dt=pi/314.0/4.0;
  Times time(ndt, dt); // ndt, dt

  /*------------------------+
  |  divergence-free field  |
  +------------------------*/
#if 1
  Comp m = Comp::u();
  for_vmijk(uvw,m,i,j,k){
    uvw[Comp::u()][i][j][k] =  uvw.zc(Comp::u(),k);
  }
  m = Comp::w();
  for_vmijk(uvw,m,i,j,k){
    uvw[Comp::w()][i][j][k] = -uvw.xc(Comp::w(),i);
  }
#else
  for(int i=0; i<uvw.ni(); i++)
    for(int j=0; j<uvw.nj(); j++)
      for(int k=0; k<uvw.nk(); k++) {
        uvw[Comp::u()][i][j][k] =  uvw.zc(Comp::u(),k);
        uvw[Comp::w()][i][j][k] = -uvw.xc(Comp::w(),i);
      }
#endif
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
#ifdef USE_VOF
  VOF T(c, f, kappa, uvw, time, solver);
#else
  CIPCSL2 T(c, f, kappa, uvw, time, solver);
  T.set_nredist(1);
  T.set_itsharpen(4);
#endif
  T.totalvol();
  boil::plot->plot(uvw, c,"uvw-c", 0);

  for(time.start(); time.end(); time.increase()) {

    T.advance();
    T.totalvol();

    if(boil::plot && time.current_step()%nint == 0) {
      boil::plot->plot(uvw,c,"uvw-c", time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-zalesak.cpp,v 1.4 2018/09/26 10:09:00 sato Exp $'/
+-----------------------------------------------------------------------------*/
