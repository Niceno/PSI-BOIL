#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#include "update_step.cpp"

#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

/* domain dimensions */
const int NX = 64;
const int NZ = 64*4;
const real LX = 0.006;
const real LZ = LX/real(NX)*real(NZ);

/* parameter for boundary or initial condition */
const real radius=1.0e-3;

/* constants */
const real gravity = 9.8;
const real pi = acos(-1.0);

/******************************************************************************/
main(int argc, char ** argv) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no() );
  Grid1D gz( Range<real>(0.0, LZ), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);           // vel
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d), step(d), sflag(d); // color function

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /* copy b.c. from p */
  press = p.shape();
  f = p.shape();

  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  g = c.shape();
  step = c.shape();
  sflag = c.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d), sapphire(d);
  vapor  .mu    (1.255e-5);
  vapor  .rho   (0.597);
  vapor  .cp    (2030*0.597);
  vapor  .lambda(0.025);
  liquid.mu    (0.28e-3);
  liquid.rho   (958.4);
  liquid.cp    (4215.9*958.4);
  liquid.lambda(0.679);
  Matter mixed(liquid, vapor, & step);
  mixed.sigma(5.9e-2);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 10000;
  const real tint = 1.0e-3;
  const real dxmin = d.dxyz_min();
  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<"\n";
  const real dt  =10.0*pow(0.5*pow(dxmin,3.0)
                 / (2.0*pi*mixed.sigma()->value()),0.5);
  const real cfl = 0.2;
  Times time(ndt, dt);
  time.print_time(false);
  time.set_coef_dec(0.75);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::ic2());

  /* color function */
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_itsharpen(8);
  conc.set_cangle(90.0);

  /* momentum equation */
  Momentum ns( uvw, xyz, time, solver, &mixed);
  //ns.convection_set(TimeScheme::forward_euler());

  /* pressure solver */
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vijk(c,i,j,k)
    c[i][j][k] = 0.0;

  const real zcent =LZ*0.8;
  for_vijk(c,i,j,k) {
    real dist=sqrt(pow(c.xc(i),2.0)+pow(c.yc(j),2.0)+pow((c.zc(k)-zcent),2.0));
    if (dist<radius*0.75) {
      c[i][j][k]=1.0;
    } else if(dist<radius*1.25) {
      int mm=8;
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
            real dist=sqrt(pow(xxc,2.0)+pow(yyc,2.0)+pow(zzc-zcent,2.0));
            if (dist<radius){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
    }
  }
  c.exchange_all();
  update_step(c, step, sflag);
  conc.init();

  boil::plot->plot(uvw,c,press,"uvw-c-press",0);

  /* set iint */
  int iint = 1;
  boil::oout<<"iint= "<<iint<<"\n";

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# DT:        " << time.dt() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "########################" << boil::endl;

    /*------------------+
    |  reset body force |
    +------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* gravity force */
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k){
      xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
    }

    /* surface tension */
    conc.tension(&xyz, mixed, step);

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* essential for moving front */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-14));

    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-6))) OMS(converged);
    p.exchange();
    ns.project(p);
    press += p;

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(d.ibody().on(i,j,k)){
        if(pmin>press[i][j][k]) pmin=press[i][j][k];
      }
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k){
      if(d.ibody().on(i,j,k)){
        press[i][j][k] -= pmin;
      } else {
        press[i][j][k] = 0.0;
      }
    }
    press.bnd_update();
    press.exchange_all();

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc.advance();
    update_step(c, step, sflag);

    /*-------------+
    |  dt control  |
    +-------------*/
    /* minimum color function */
    time.control_dt(ns.cfl_max(), cfl, dt);

    /*--------------+
    |  output data  |
    +--------------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,c,press,"uvw-c-press",iint);
      iint = int(time.current_time()/tint) + 1;
    }

  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}
