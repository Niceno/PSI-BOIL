#include "Include/psi-boil.h"
#include <fstream>

/* domain dimensions (given by problem) */
//const real LX =   0.5;
const real LX =   0.25;
const real LY =   0.016;
const real LZ =   0.008;

/* computed parameters */
//const int NX = 4;
const int NX = 125;
const int NY = 4;
const int NZ = 4;

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
  //Grid1D gx( Range<real>( 0,LX), NX, Periodic::no() );
  Grid1D gx( Range<real>( 0,LX), NX, Periodic::yes() );
  Grid1D gy( Range<real>( 0,LY), NY, Periodic::yes());
  Grid1D gz( Range<real>( 0,LZ), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 500;
  //const int ndt = 1;
  const int nint = 50;
  const real dt  = 0.25 * LX / real(NX);
  Times time(ndt, dt); 
	
  OPR(  NX );
  OPR(  dt );
  OPR( ndt );

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // vel
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    //uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet(), 1.0,0.0,0.0 ) );
    //uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet()   ) );
    //uvw.bc(m).add( BndCnd( Dir::imax(), BndType::inlet(),-1.0,0.0,0.0 ) );
    //uvw.bc(m).add( BndCnd( Dir::imin(), BndType::outlet()   ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }
  
  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  //c.bc().add( BndCnd( Dir::imin(), BndType::inlet(),1.0 ) );
  //c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  //c.bc().add( BndCnd( Dir::imax(), BndType::inlet(),1.0 ) );
  //c.bc().add( BndCnd( Dir::imin(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  Comp m=Comp::u();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=1.0;

  m=Comp::v();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=0.0;

  m=Comp::w();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=0.0;

  uvw.exchange_all();

  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

#if 1
  for_vijk(c,i,j,k) {
    if( 0.10 < c.xc(i) && c.xc(i) < 0.20 ){
          c[i][j][k]=1.0;
    }
  }
#else
  for_vijk(c,i,j,k) {
    c[i][j][k]=std::min(1.0,c.xc(i)*5.0);
    //c[i][j][k]=sin(c.xc(i)*5.0);
  }
#endif
  
  c.exchange();
  boil::plot->plot(uvw,c, "uvw-c-init0", 0);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (1.0000e-5);
  air  .rho   (1.1768e+0);
  water.mu    (1.0000e-3);
  water.rho   (1.0000e+3);

  //Matter mixed(air, water, c);
  Matter mixed(water, air, c);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  CIPCSL2 conc  (c, g, kappa, uvw, time, solver);
#if 1
  std::ofstream fout0;
  fout0.open("init-profile.txt");
  for_vi(c,i) {
       fout0 << c.xc(i) << "  " << c[i][1][1] << "\n";
  }
  fout0.close();
#endif

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "########################" << boil::endl;
	  
    /*---------------------------+
    |  fully explicit with conc  |
    +---------------------------*/
    conc.new_time_step();
    //conc.convection();
    conc.advance();
    //conc.sharpen();
    conc.totalvol();

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
    }

  }

#if 1
  std::ofstream fout;
  fout.open("profile.txt");
  for_vi(c,i) {
       fout << c.xc(i) << "  " << g[i][1][1] << "\n";
  }
  fout.close();
  return 0;
#endif

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-CIPCSL2-1d.cpp,v 1.3 2018/09/26 10:06:18 sato Exp $'/
+-----------------------------------------------------------------------------*/
