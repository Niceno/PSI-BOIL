#include "Include/psi-boil.h"
#include <fstream>
#define USE_VOF
//#define GHOST

#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

/* computed parameters */
//const int NX = 50;
const int gLevel = 4;
const int NX = 32*gLevel;
const int NY = 32*gLevel;
const int NZ = 32*gLevel;
//const int NZ = 50;

/* domain dimensions (given by problem) */
const real LX =   1.0;
const real LY =   1.0;
const real LZ =   1.0;
const real radius = LX/4.0;


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
#if 0
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no() );
  Grid1D gy( Range<real>(-0.5*LY,0.5*LY), NY, Periodic::no() );
  Grid1D gz( Range<real>(-0.5*LZ,0.5*LZ), NZ, Periodic::no() );
#else
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no() );
  Grid1D gy( Range<real>(-0.5*LY,0.5*LY), NY, Periodic::yes() );
  Grid1D gz( Range<real>(-0.5*LZ,0.5*LZ), NZ, Periodic::no() );
#endif
  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // vel
  Vector xyz(d); // body force
  Scalar c  (d), g  (d), kappa(d); // concentration
  Scalar press  (d), p(d), f  (d);


  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
#if 0
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
#else
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
#endif
  }

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );
#if 0
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
#else
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
#endif

  press.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  Matter air(d), water(d);
  air  .mu    (2.0000e-5);
  air  .rho   (1.2500e+0);
  water.mu    (1.0000e-3);
  water.rho   (1.0000e+3);

  Matter mixed(water, air, &c);
  mixed.sigma (0.07);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  real dxmin = d.dxyz_min();
  const real dt = 1.0 * pow(air.rho()->value()*pow(dxmin,3.0)
                           /(2.0*3.1415*mixed.sigma()->value()),0.5);
#if 1
  const real tend = 1.0;
  const int ndt = tend/dt;
  const int nint = ndt;
#else
  const int  ndt  = 351;
  const int nint  = 50;
  const int nint2 = 350;
#endif

  Times time(ndt, dt); 

  
  bool restart = false;
  int ts=0;
  std::fstream input;
  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    restart=true;
  }

  if( restart ) {
    real t,dtf;
    input >> ts;
    input >> t;
    time.first_step(ts);
    time.current_time(t);
    uvw.load("uvw",ts);
    press.load("press",ts);
    c.load("c",ts);
  } else {
    /*--------------------+
    |  initial condition  |
    +--------------------*/
    Comp m=Comp::u();
    for_vmijk(uvw,m,i,j,k)
      uvw[m][i][j][k]=0.0;
  
    m=Comp::v();
    for_vmijk(uvw,m,i,j,k)
      uvw[m][i][j][k]=0.0;

    m=Comp::w();
    for_vmijk(uvw,m,i,j,k)
      uvw[m][i][j][k]=0.0;
    uvw.exchange();

    for_vijk(c,i,j,k) 
      c[i][j][k] = 0.0;
  
    const real xcent = 0.0;
    const real ycent = 0.0;
    const real zcent = 0.0;
    for_vijk(c,i,j,k) {
      real dist=sqrt(pow(c.xc(i)-xcent,2.0)
                    +pow(c.yc(j)-ycent,2.0)+pow(c.zc(k)-zcent,2.0));
      //real dist=sqrt(pow(c.xc(i)-xcent,2.0)
      //              +pow(c.zc(k)-zcent,2.0));
      if (dist<radius*0.80) {
        c[i][j][k]=1.0;
      } else if(dist<radius*1.2) {
        int mm=20;
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
              //real dist=sqrt(pow(xxc-xcent,2.0)
              //              +pow(zzc-zcent,2.0));
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
    boil::plot->plot(uvw,c,press, "uvw-c-press", 0);
  }

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());
  Momentum ns( uvw, xyz, time, solver, &mixed);
  Pressure pr(press, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );

#ifdef USE_VOF
  VOF conc  (c, g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c, g, kappa, uvw, time, solver);
#endif
  conc.totalvol();


  for(time.start(); time.end(); time.increase()) {

    //reset
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* surface tension */
#ifndef GHOST
    conc.tension(&xyz, mixed);
#endif

    /* essential for moving front */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();
#ifndef GHOST
    ns.grad(press);
#endif
    ns.solve(ResRat(1e-8));

#ifdef GHOST
    f = 0.0;
    conc.curvature();
    pr.ghost(c,kappa);
#else
    p=0.0;
#endif

    //if (multigrid.vcycle(ResRat(1e-4))) OMS(converged);
    pr.solve(ResRat(1e-8));

#ifdef GHOST
    ns.project_ghost(press,c,kappa);
#else
    ns.project(p);
    press += p;
#endif

    conc.advance();
    conc.totalvol();

    /*---------------+
    |  post process  |
    +-------------- */
    // pressure
    real dltp = -1.0e+300;
    for_vi(press,i) {
      if(approx(press.xc(i),gx.xc(gx.ncell()),1.0e-6)) {
        for_vj(press,j) {
          if(approx(press.yc(j),gy.xc(gy.ncell()),1.0e-6)) {
            for_vk(press,k) {
              if(approx(press.zc(k),gz.xc(gz.ncell()),1.0e-6)) {
                dltp = press[i][j][k];
              }
            }
          }
        }
      }
    }
    boil::cart.max_real(&dltp);

    for_vijk(press,i,j,k)
      press[i][j][k] -= dltp;

    int icount=0;
    real psum=0.0;
    real perr=0.0;
    const real pexact=mixed.sigma()->value()/radius;
    for_vijk(press,i,j,k) {
      if(c[i][j][k]>=0.99){
        psum += press[i][j][k];
        icount ++;
        perr=pow((press[i][j][k]-pexact),2.0);
      }
    }
    boil::cart.sum_real(&psum);
    boil::cart.sum_real(&perr);
    boil::cart.sum_int(&icount);

    if (icount==0) {
      boil::oout<<"Error: no liquid cells\n";
      exit(0);
    } else {
      perr=perr/(real(icount)*pow(pexact,2.0));
      perr=sqrt(perr);
      boil::oout<<"pressure= "<<time.current_time()<<" "<<psum/real(icount)<<" "
               <<pexact<<" "<<perr<<" "<<icount<<"\n";
    }

    //velocity
    real umax=0.0, vmax=0.0, wmax=0.0;
    real usum=0.0, vsum=0.0, wsum=0.0;
    icount=0;
    for_vijk(c,i,j,k) {
      real utmp = 0.5*(uvw[Comp::u()][i][j][k]+uvw[Comp::u()][i+1][j][k]);
      real vtmp = 0.5*(uvw[Comp::v()][i][j][k]+uvw[Comp::v()][i][j+1][k]);
      real wtmp = 0.5*(uvw[Comp::w()][i][j][k]+uvw[Comp::w()][i][j][k+1]);
      if (umax<fabs(utmp)) umax=fabs(utmp);
      if (vmax<fabs(vtmp)) vmax=fabs(vtmp);
      if (wmax<fabs(wtmp)) wmax=fabs(wtmp);
      usum += utmp;
      vsum += vtmp;
      wsum += wtmp;
      icount++;
    }
    boil::cart.sum_real(&usum);
    boil::cart.sum_real(&vsum);
    boil::cart.sum_real(&wsum);
    boil::cart.sum_int(&icount);
    usum /= real(icount);
    vsum /= real(icount);
    wsum /= real(icount);
    boil::cart.max_real(&umax);
    boil::cart.max_real(&vmax);
    boil::cart.max_real(&wmax);
    boil::oout<<"velocity= "<<time.current_time()<<" "
             <<umax<<" "<<vmax<<" "<<wmax<<" "
             <<usum<<" "<<vsum<<" "<<wsum<<"\n";

    // kappa
    real kappa_min= 1.0e+16;
    real kappa_max=-1.0e+16;
    for_vijk(kappa,i,j,k) {
      if (kappa_min>kappa[i][j][k]) kappa_min=kappa[i][j][k];
      if (kappa_max<kappa[i][j][k]) kappa_max=kappa[i][j][k];
    }
    boil::cart.min_real(&kappa_min);
    boil::cart.max_real(&kappa_max);
    boil::oout<<"kappa_min_max= "<<kappa_min<<" "<<kappa_max<<"\n";


    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, press,"uvw-c-press",  time.current_step());
      boil::plot->plot(xyz, c, kappa,"xyz-c-kappa",  time.current_step());
    }
#if 0
    if(time.current_step() % nint2 == 0) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      c    .save("c",     time.current_step());
      std::fstream output;
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output.close();
    }
#endif
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-CIPCSL2-1d.cpp,v 1.3 2018/09/26 10:06:18 sato Exp $'/
+-----------------------------------------------------------------------------*/
