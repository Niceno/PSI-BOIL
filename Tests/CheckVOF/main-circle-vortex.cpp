#include "Include/psi-boil.h"
#define USE_VOF

/* domain dimensions (given by problem) */
const real LX =   1.0;
const real LZ =   0.05; //0.004;

/* computed parameters */
const int NX = 128;
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
  Grid1D gx( Range<real>( 0.0,LX), NX, Periodic::no() );
  Grid1D gz( Range<real>( 0.0,LZ), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );

  /* copy b.c. from p */
  
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );

  g = c.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (5.0000e-4);
  air  .rho   (5.0000e+2);
  water.mu    (1.0000e-3);
  water.rho   (1.0000e+3);

  Matter mixed(water, air, &c);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 8000;
  const int  nint = 1000;
  const real dt  = 0.001;

  Times time(ndt, dt);

  OPR(  NX );
  OPR(  dt );
  OPR( ndt );

  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

  const real radius=0.15;
  const real xcent=0.5;
  const real ycent=0.75;
  for_vijk(c,i,j,k) {
    real dist=pow(c.xc(i)-xcent,2.0)+pow(c.yc(j)-ycent,2.0);
    if (dist<pow(radius*0.75,2.0)) {
      c[i][j][k]=1.0;
    } else if(dist<pow(radius*1.25,2.0)) {
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
            real dist=pow(xxc-xcent,2.0)+pow(yyc-ycent,2.0);
            if (dist<pow(radius,2.0)){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
    }
  }

  
  c.exchange();
  boil::plot->plot(uvw, c, "uvw-c",  0);

#ifdef USE_VOF
  VOF conc  (c,  g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
#endif
  conc.totalvol();

  for(time.start(); time.end(); time.increase()) {

    real pi=acos(-1.0);
    real period=8.0;
    real t=time.current_time();
    Comp m = Comp::u();
    for_avmijk(uvw,m,i,j,k){
      real x=uvw.xc(m,i);
      real y=uvw.yc(m,j);
      uvw[m][i][j][k] = -sin(pi*x)*sin(pi*x)
                       *sin(2.0*pi*y)*cos(pi*t/period);
    }
    m = Comp::v();
    for_avmijk(uvw,m,i,j,k){
      real x=uvw.xc(m,i);
      real y=uvw.yc(m,j);
      uvw[m][i][j][k] =  sin(2*pi*x)
                         *sin(pi*y)*sin(pi*y)*cos(pi*t/period);
    }
    uvw.exchange_all();

    // check volume convergence
    for_vijk(c,i,j,k) {
      real vol_conv = (-uvw[Comp::u()][i][j][k]+uvw[Comp::u()][i+1][j][k]
                       -uvw[Comp::v()][i][j][k]+uvw[Comp::v()][i][j+1][k]
                       -uvw[Comp::w()][i][j][k]+uvw[Comp::w()][i][j][k+1])
                      *time.dt()/c.dV(i,j,k);
      if (fabs(vol_conv)>1e-12) std::cout<<"vol_conv= "<<vol_conv<<" "
                                        <<i<<" "<<j<<" "<<k<<"\n";
    }

    conc.advance();
    conc.totalvol();

    if(time.current_step() % (nint)==0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-circle-vortex.cpp,v 1.4 2018/09/26 10:09:00 sato Exp $'/
+-----------------------------------------------------------------------------*/
