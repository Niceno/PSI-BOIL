/*----------------------------------------------------------------------------+
|  Rising 2D bubble in stagnant liquid                                        |
|  https://doi.org/10.1002/fld.1934                                           |
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"
#define USE_VOF   // activate this line to use VOF method

#ifndef USE_VOF
void update_step(const Scalar & c, Scalar & step, Scalar & sflag);
#endif

const int gLevel = 4;  // h = 1/160
const int NX = 40*gLevel;
const int NY = 2*NX;
const int NZ = 4;

const real LX = 1.0;
const real LY = 2.0;
const real LZ = LX/real(NX)*real(NZ);
const real radius = 0.25;  // diameter = 0.5 m

const real gravity = -0.98; //[m/s2]

/******************************************************************************/
int main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no(),BndGrid::symmetry(),BndGrid::symmetry());
  Grid1D gy( Range<real>(0, LY), NY, Periodic::no() );
  Grid1D gz( Range<real>(-0.5*LZ,0.5*LZ), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::ic2());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d); // concentration
#ifndef USE_VOF
  Scalar step(d), sflag(d);
#endif

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::symmetry() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  press = p.shape();
  f = p.shape();
  g = p.shape();
  kappa = p.shape();
  
  c.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
#ifndef USE_VOF
  step = c.shape();
  sflag = c.shape();
#endif

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), oil(d);
  air.rho   (1.0);
  air.mu    (0.1);
  oil.rho   (1000.0);
  oil.mu    (10.0); //(d)
#ifdef USE_VOF
  Matter mixed(oil, air, & c);
#else
  Matter mixed(oil, air, & step);
#endif
  mixed.sigma(1.96);

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 10.0 * pow(air.rho()->value()*pow(dxmin,3.0)
                        /(2.0*boil::pi*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 100000;
  Times time(ndt, dt);
  const real tint = 0.1;
  const int nint= 10000;
  const real cfl_limit=0.25;

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.max_cycles(10);
  multigrid.min_cycles(4);

  /*--------------------+
  |  initial condition  |
  +--------------------*/
#if 1
  for_vijk(c,i,j,k)
    c[i][j][k] = 1.0;

  const real ycent =0.5;
  for_vijk(c,i,j,k) {
    real dist=sqrt(pow(c.xc(i),2.0)+pow(c.yc(j)-ycent,2.0)+pow(c.zc(k),2.0));
    if (dist<radius*0.9) {
      c[i][j][k]=0.0;
    } else if(dist<radius*1.1) {
      int mm=16;
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
            real dist=sqrt(pow(xxc,2.0)+pow(yyc-ycent,2.0)+pow(zzc,2.0));
            if (dist>radius){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
    }
  }
  c.exchange_all();

  /* output initial condition */
  boil::plot->plot(uvw,c,press, "uvw-c-press",0);

#else
  int nrestart=3500;
  uvw.load("uvw", nrestart);
  press.load("press", nrestart);
  c.  load("c",   nrestart);
  uvw.exchange_all();
  p.exchange_all();
  c.exchange_all();
#endif

#ifdef USE_VOF
  VOF conc (c,   g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc (c,   g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
  conc.front_minmax();
#endif

#ifndef USE_VOF
  update_step(c, step, sflag);
#endif

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  Time loop  |
  +------------*/

  for(time.start(); time.end(); time.increase()) {

    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* add gravity force */ 
    Comp m = Comp::v();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }
    xyz.exchange();

    /* add surface tension force */
#ifdef USE_VOF
    conc.tension(&xyz, mixed);
#else
    conc.tension(&xyz, mixed, step);
#endif

    /* solve momentum equation */
    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(0.001));

    p = 0.0;
    multigrid.vcycle(ResRat(1e-3));
    ns.project(p);
    press += p;
    press.exchange();

    /* advance VOF or color function */
    conc.new_time_step();
    conc.advance();
#ifndef USE_VOF
    update_step(c, step, sflag);
#endif
    conc.totalvol();
    conc.front_minmax();

    /* dt control */
    time.control_dt(ns.cfl_max(), cfl_limit, dt);

    /* output field data */
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,c,press, "uvw-c-press",iint);

      // subtract gravity just for visualization
      m = Comp::v();
      for_vmijk(xyz,m,i,j,k) {
        xyz[m][i][j][k] -= gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
      }
      xyz.exchange();
      boil::plot->plot(xyz,c,kappa, "xyz-c-kappa",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /* backup field data for restart */
    if( time.current_step() % nint == 0 ) {
      uvw.save("uvw", time.current_step());
      press.save("press", time.current_step());
      c.save("c", time.current_step());
    }

    /*---------+
    |  exit ?  |
    +---------*/
    if(conc.get_ymaxft()>=LY*0.95){
       std::cout<<"Bubble reaches to the top boundary. Exit.";
       exit(0);
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
	
/*-----------------------------------------------------------------------------+
 '$Id: main-bubble-3d.cpp,v 1.5 2009/07/01 14:18:53 niceno Exp $'/
+-----------------------------------------------------------------------------*/
