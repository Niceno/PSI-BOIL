/*----------------------------------------------------------------------------+
|  Rising air bubble in stagnant oil                                          |
|  Publication: Y.Sato, B.Niceno, Int J Numer Meth Fluids, 70 (2012) 441-467  |
|  Case (d)-de/20 in Fig.23                                                   |
+----------------------------------------------------------------------------*/
#include "Include/psi-boil.h"
#define USE_VOF

void update_step(const Scalar & c, Scalar & step, Scalar & sflag);

const int NX= 120;
const int NZ= NX*2;
const real RB = 0.5 * 0.0261;  // diameter = 0.0261 m
const real LX = RB*12.0;
const real LZ = LX*NZ/NX; 

const real gravity = -9.8;

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
  Matter air(d), oil(d);
  air  .rho   (1.205);
  air.mu    (1.82e-5);
  oil.rho   (1350.);
  oil.mu    (0.7745); //(d)

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
#ifndef USE_VOF
  step = c.shape();
  sflag = c.shape();
#endif

#ifdef USE_VOF
  Matter mixed(oil, air, & c);
#else
  Matter mixed(oil, air, & step);
#endif
  mixed.sigma(0.078);

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = d.dxyz_min();
  const real dt  = 5.0 * pow(air.rho()->value()*pow(dxmin,3.0)
                        /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 10000;
  Times time(ndt, dt);
  const real tint = 0.1;
  const int nint= 24000;
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

  const real radius=RB;
  const real zcent =LZ*0.1;
  for_vijk(c,i,j,k) {
    real dist=sqrt(pow(c.xc(i),2.0)+pow(c.yc(j),2.0)+pow((c.zc(k)-zcent),2.0));
    if (dist<radius*0.75) {
      c[i][j][k]=0.0;
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
#endif

  //boil::plot->plot(uvw,c, press, "uvw-c-press",0);
  boil::plot->plot(uvw,c,press, "uvw-c-press",0);
  conc.front_minmax();

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
 
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] = gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
    }
    xyz.exchange();


#ifdef USE_VOF
    conc.tension(&xyz, mixed);
#else
    conc.tension(&xyz, mixed, step);
#endif

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.grad(press);
    ns.convection();
    ns.solve(ResRat(0.001));

    p = 0.0;
    //multigrid.max_cycles(10);
    multigrid.vcycle(ResRat(1e-3));
    ns.project(p);
    press += p;
    press.exchange();

    /* advance */
    conc.new_time_step();
    conc.advance();
#ifndef USE_VOF
    update_step(c, step, sflag);
#endif
    conc.totalvol();
    conc.front_minmax();

    /* dt control */
    time.control_dt(ns.cfl_max(), cfl_limit, dt);

    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,c,press, "uvw-c-press",iint);

      // subtract gravity just for visualization
      m = Comp::w();
      for_vmijk(xyz,m,i,j,k) {
        xyz[m][i][j][k] -= gravity*xyz.dV(m,i,j,k)*mixed.rho(m,i,j,k);
      }
      xyz.exchange();
      boil::plot->plot(xyz,c,kappa, "xyz-c-kappa",iint);
      iint = int(time.current_time()/tint) + 1;
    }
    if( time.current_step() % nint == 0 ) {
      uvw.save("uvw", time.current_step());
      press.  save("press",   time.current_step());
      c.  save("c",   time.current_step());
    }

    /*---------+
    |  exit ?  |
    +---------*/
    if(conc.get_zmaxft()>=LZ*0.95){
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
