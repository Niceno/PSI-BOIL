/*----------------------+
| Rising bubble         |
+----------------------*/
#include "Include/psi-boil.h"

#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

void update_step(const Scalar & c, Scalar & step, Scalar & sflag);

const int Level=1;  // =1,2,4
const int NX= 32*Level;
const int NZ= NX*2;
const int nint = 50*Level;  //plot every nint time steps
const int ndt = 2500;   //run for ndt time steps

const real RB = 0.5 * 3.0e-3;
const real LX = RB*8.0;
const real LZ = LX*NZ/NX; 

const real gravity=9.8;

/* make spherical bubbles */
void dispersed_sphere_init(real, real, real, real, Scalar &, real,
                           const real);

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
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 
             Range<real>( LX/NX,  LX/NX ),
              NX, Periodic::no());

  Grid1D gz( Range<real>(0.0, LZ), 
             Range<real>( LZ/NZ,  LZ/NZ ),
              NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .rho   (1.205);
  air.mu    (1.82e-5);
  water.rho  (998.2);
  water.mu   (1.0e-3);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d), step(d), sflag(d); // concentration
  Scalar kappa(d);
  Scalar press(d);
  Scalar rgnid(d);  //to track regions

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
  press = p.shape();
  
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );
  step = c.shape();
  kappa = c.shape();
  sflag = c.shape();
  rgnid = c.shape();

  Matter mixed(water, air, &step);
  mixed.sigma(0.072);
  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = std::min(LX/NX,LZ/NZ);
  const real dt  = 5.0 * pow(0.5*pow(dxmin,3.0)
                           /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
//  const int ndt = 80000;
  Times time(ndt, dt);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());

  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.max_cycles(30);
  multigrid.min_cycles(3);
  multigrid.stop_if_diverging(false);

  /*-----------------------------------------------------------------------------+
  |  Tracking algorithm object
  +-----------------------------------------------------------------------------*/
  Floodfill idrgns(c, rgnid, &uvw, time);
  idrgns.set_out_freq(100);
  boil::oout<<"main:Floodfill:output frequency= "<<idrgns.get_out_freq()<<"\n";

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vijk(c,i,j,k) {
    c[i][j][k] = 0.0;  // c = 0.0 air    1.0 water
    if(c.zc(k) < 0.80*LZ) c[i][j][k] = 1.0; 
  }
  dispersed_sphere_init (LX* 0.0, LX*0.0, LZ*0.09, 1.0, c , 0.0, RB);
  dispersed_sphere_init (LX* 0.0, LX*0.0, LZ*0.32, 2.0, c , 0.0, RB);  
  dispersed_sphere_init (LX* 0.0, LX*0.0, LZ*0.54, 0.7, c , 0.0, RB);
  c.exchange_all();

  boil::plot->plot(c, "c",0);
//  exit(1);

  //ColorFunction  conc  (c,   g, uvw, time, solver); 
  CIPCSL2 conc (c,   g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
  conc.front_minmax();
  update_step(c, step, sflag);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;
 
    /* advance */
    conc.new_time_step();
    conc.advance();
    conc.totalvol();
    conc.front_minmax();
    conc.tension(&xyz, mixed, step);
    update_step(c, step, sflag);

    const Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] -= gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.grad(press);
    ns.convection();
    ns.solve(ResRat(1e-4));

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-3));
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();
    
    /* Identify region with tracking algorithm */
    idrgns.identify_regions();

    /* dt control */
    time.control_dt(ns.cfl_max(), 0.25, dt);

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, rgnid, "uvw-c-rgnid",time.current_step()/nint);
    }
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}

/******************************************************************************/
void update_step(const Scalar & c, Scalar & step, Scalar & sflag){
  const real phisurf=0.5;
  for_avijk(sflag,i,j,k) {
    sflag[i][j][k]=0;
  }
  for_vijk(c,i,j,k) {
    if(c[i][j][k]>=phisurf)
      sflag[i][j][k]=1;
  }
  /* i-direction */
  for(int i=c.si()-1; i<=c.ei(); i++){
    for_vjk(c,j,k){
       if((c[i][j][k]-phisurf)*(c[i+1][j][k]-phisurf)<=0.0){
          sflag[i  ][j][k]=2;
          sflag[i+1][j][k]=2;
       }
    }
  }
  /* j-direction */
  for(int j=c.sj()-1; j<=c.ej(); j++){
    for_vik(c,i,k){
      if((c[i][j][k]-phisurf)*(c[i][j+1][k]-phisurf)<=0.0){
          sflag[i][j  ][k]=2;
          sflag[i][j+1][k]=2;
       }
    }
  }
  /* k-direction */
  for(int k=c.sk()-1; k<=c.ek(); k++){
    for_vij(c,i,j){
       if((c[i][j][k]-phisurf)*(c[i][j][k+1]-phisurf)<=0.0){
          sflag[i][j][k  ]=2;
          sflag[i][j][k+1]=2;
       }
    }
  }
  sflag.exchange_all();
  for_avijk(c,i,j,k){
    if(sflag[i][j][k]==2){
      step[i][j][k]=c[i][j][k];
    } else {
      if(c[i][j][k]<0.5){
        step[i][j][k]=0.0;
      } else {
        step[i][j][k]=1.0;
      }
    }
  }
}

/*****************************************************************************/
/* function to make drops or bubble */
void dispersed_sphere_init (real xcent, real ycent, real zcent, real rad_fact, 
     Scalar & c , real cf, const real radius ) {
  for_vijk(c,i,j,k) {
    real dist=pow((c.xc(i)-xcent),2.0)+pow((c.yc(j)-ycent),2.0)+pow((c.zc(k)-zcent),2.0);
    real rad = radius * rad_fact;
    if (dist<pow(rad*0.75,2.0)) {
      c[i][j][k]=cf;
    } 
    else if(dist<pow(rad*1.25,2.0)) {
      int mm=8;
      real x0=c.xn(i);
      real y0=c.yn(j);
      real z0=c.zn(k);
      real ddx=c.dxc(i)/real(mm);
      real ddy=c.dyc(j)/real(mm);
      real ddz=c.dzc(k)/real(mm);
      int itmp=0;
      for (int ii=0; ii<mm; ii++){
        for (int jj=0; jj<mm; jj++){
          for (int kk=0; kk<mm; kk++){
            real xxc=x0+0.5*ddx+real(ii)*ddx;
            real yyc=y0+0.5*ddy+real(jj)*ddy;
            real zzc=z0+0.5*ddz+real(kk)*ddz;

            real dist=pow(xxc-xcent,2.0)+pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0);
            if (dist<pow(rad,2.0)){  //
                itmp=itmp+1;
            }
          }
        }
      }
      if (fabs(c[i][j][k]-cf) > 1.0e-8)
        c[i][j][k]=fabs(c[i][j][k]-real(itmp)/real(mm*mm*mm));
    }
  }
  c.exchange_all();
}

/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.1 2018/02/16 19:06:00 sato Exp $'/
+-----------------------------------------------------------------------------*/
