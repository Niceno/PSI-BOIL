#include "Include/psi-boil.h"
#include <fstream>
#define USE_VOF

#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}


/* domain dimensions (given by problem) */
//const real LX =   0.5;
const real LX =   0.2;
const real LY =   0.2;
const real LZ =   0.2;

/* computed parameters */
//const int NX = 4;
const int NX = 100;
const int NY = 100;
const int NZ = 100;

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
  Grid1D gx( Range<real>( 0,LX), NX, Periodic::yes() );
//  Grid1D gx( Range<real>( 0,LX), Range<real>(LX/(2.0 * NX),LX/(2.0 * NX)), NX, Periodic::yes() );
  Grid1D gy( Range<real>( 0,LY), NY, Periodic::yes());
//  Grid1D gy( Range<real>( 0,LY), Range<real>(LY/(2.0 * NY),LY/(2.0 * NY)), NY, Periodic::yes() );
  Grid1D gz( Range<real>( 0,LZ), NZ, Periodic::yes() );
//  Grid1D gz( Range<real>( 0,LZ), Range<real>(LZ/(2.0 * NZ),LZ/(2.0 * NZ)), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 4000;
  //const int ndt = 1;
  const int nint = 50;
  const real dt  = 0.25 * LX / real(NX);
  Times time(ndt, dt); 
	
  //OPR(  NX );
  //OPR(  dt );
  //OPR( ndt );

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // vel
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
#if 0
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
#endif
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
    uvw[m][i][j][k]=1.0;

  m=Comp::w();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=1.0;

  uvw.exchange_all();

  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

#if 1
  for_vijk(c,i,j,k) {
    if( 0.02 < c.xc(i) && c.xc(i) < 0.06 && 0.02 < c.yc(j) && c.yc(j) < 0.06 && 0.02 < c.zc(k) && c.zc(k) < 0.06){
          c[i][j][k]=1.0;
    }

  }

  //std::cout<<"c50 "<<c[50][1][1]<<" "<<"c51 "<<c[51][1][1]<<"\n"; 

  //std::cout<<"c75 "<<c[75][1][1]<<" "<<"c76 "<<c[76][1][1]<<"\n"; 

#else
  for_vijk(c,i,j,k) {
    c[i][j][k]=std::min(1.0,c.xc(i)*5.0);
    //c[i][j][k]=sin(c.xc(i)*5.0);
  }
#endif
  
  c.exchange_all();
  c.bnd_update();
  boil::plot->plot(uvw,c, "uvw-c-init0", 0);


  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

#ifdef USE_VOF
  VOF conc  (c, g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c, g, kappa, uvw, time, solver);
#endif

#if 0
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
   // new_time_step();

    //std::cout<<"c50 before advance "<<c[50][1][1]<<" "<<"c51 before advance "<<c[51][1][1]<<"\n";
    //std::cout<<"c75 before advance "<<c[75][1][1]<<" "<<"c76 before advance "<<c[76][1][1]<<"\n"; 
    conc.advance();
    conc.totalvol();
    //std::cout<<"c50 after advance "<<c[50][1][1]<<" "<<"c51 after advance "<<c[51][1][1]<<"\n"; 
    //std::cout<<"c75 after advance "<<c[75][1][1]<<" "<<"c76 after advance "<<c[76][1][1]<<"\n"; 
    //conc.convection();
    //conc.advance();
    //conc.sharpen();
    //conc.totalvol();

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
    }

  }

     
#if 0
  std::ofstream fout;
  fout.open("profile.txt");
  for_vi(c,i) {
       fout << c.xc(i) << "  " << c[i][1][1] << "\n";
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
