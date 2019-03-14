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
  Grid1D gy( Range<real>( 0,LY), NY, Periodic::no() );
  Grid1D gz( Range<real>( 0,LZ), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 400;
  const int nint = 100;
  const real dt  = 0.25 * LX / real(NX);
  Times time(ndt, dt); 
	
  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // vel
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
#if 1
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
  }
#endif
  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  Comp m=Comp::u();
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=1.0;

  m=Comp::v();
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=0.0;

  m=Comp::w();
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=0.0;
  uvw.exchange();

  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

  const real radius = LZ/4.0 - 0.5*LZ/real(NZ);
  const real xcent = LX/2.0;
  const real ycent = LX/2.0;
  const real zcent = LZ/4.0;
  for_vijk(c,i,j,k) {
    real dist=sqrt(pow(c.xc(i)-xcent,2.0)
                  +pow(c.yc(j)-ycent,2.0)+pow(c.zc(k)-zcent,2.0));
    if (dist<radius*0.80) {
      c[i][j][k]=1.0;
    } else if(dist<radius*1.2) {
      int mm=10;
      real x0=d.xn(i);
      real y0=d.yn(j);
      real z0=d.zn(k);
      real ddx=d.dxc(i)/real(mm);
      real ddy=d.dyc(j)/real(mm);
      real ddz=d.dzc(k)/real(mm);
      int itmp=0;
      for (int ii=0; ii<mm; ii++){
        for (int jj=0; jj<mm; jj++){
          for (int kk=0; kk<mm; kk++){
            real xxc=x0+0.5*ddx+real(ii)*ddx;
            real yyc=y0+0.5*ddy+real(jj)*ddy;
            real zzc=z0+0.5*ddz+real(kk)*ddz;
            real dist=sqrt(pow(xxc-xcent,2.0)
                          +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
            if (dist<radius){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
    }
  }
  c.bnd_update();
  c.exchange_all();
  boil::plot->plot(uvw,c, "uvw-c", 0);


  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

#ifdef USE_VOF
  VOF conc  (c, g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c, g, kappa, uvw, time, solver);
#endif
  conc.totalvol();


  for(time.start(); time.end(); time.increase()) {

    conc.advance();
    conc.totalvol();

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
    }

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-CIPCSL2-1d.cpp,v 1.3 2018/09/26 10:06:18 sato Exp $'/
+-----------------------------------------------------------------------------*/
