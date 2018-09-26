#include "Include/psi-boil.h"
#define GHOST

/* computed parameters */
const int NX = 30;
const int NY = 4;

/* domain dimensions (given by problem) */
const real LX =   0.06;
const real LY =   LX/real(NX)*real(NY);

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
  //real zcent = 0.48*LX;
  real zcent = 0.3*LX;
  Grid1D gx( Range<real>( -0.5*LX,0.5*LX), NX, Periodic::no() );
  Grid1D gz( Range<real>( -0.5*LX+zcent,0.5*LX+zcent), NX, Periodic::no() );
  Grid1D gy( Range<real>( 0.0,LY), NY, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);
  //Body floor("floor.stl");
  //Domain d(gx, gy, gz, & floor);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar press  (d), f  (d);
#ifndef GHOST
  Scalar p (d);
#endif
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }
  
  press.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );

  g = c.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (5.0000e-4);
  air  .rho   (5.0000e+2);
  water.mu    (1.0000e-3);
  water.rho   (1.0000e+3);

  Matter mixed(water, air, & c);
  mixed.sigma(2.3610e-2);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 10000;
  const int  nint = 1000;
  const real dxmin = std::min(LX/NX,LY/NY);
  const real dt  = 0.1*pow(air.rho()->value()*pow(dxmin,3.0)/(2.0*3.1415*mixed.sigma()->value()),0.5);

  Times time(ndt, dt);

  OPR(  NX );
  OPR(  dt );
  OPR( ndt );


  
  for_vijk(c,i,j,k) {
    c[i][j][k] = 1.0;  // bubble
    //c[i][j][k] = 0.0;  // droplet
  }

  const real radius=0.02;
  const real xcent=0.0;
  const real ycent=-0.00;
  for_vijk(c,i,j,k) {
    real dist=pow(c.xc(i)-xcent,2.0)+pow(c.zc(k)-zcent,2.0);
    if (dist<pow(radius*0.75,2.0)) {
      c[i][j][k]=0.0;  // bubble
      //c[i][j][k]=1.0;  // droplet
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
            real zzc=z0+0.5*ddz+real(kk)*ddz;
            real dist=pow(xxc-xcent,2.0)+pow(zzc-zcent,2.0);
            if (dist<pow(radius,2.0)){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=1.0-real(itmp)/real(mm*mm*mm);  // bubble
      //c[i][j][k]=real(itmp)/real(mm*mm*mm);  // droplet
    }
  }

  
  c.exchange();
  boil::plot->plot(uvw, c, press, "uvw-c-press",  0);

  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.front_minmax();
  conc.totalvol();
  conc.set_cangle(0.0);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);

  Pressure pr(press, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);


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
    conc.advance();
    conc.totalvol();

    /* essential for moving front */
    ns.discretize();
    pr.discretize();

    /* momentum */
    ns.new_time_step();

    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* surface tension */
#ifndef GHOST
    conc.tension(&xyz, mixed);
    std::cout<<"cangle= "<<conc.get_cangle()<<"\n";
    //boil::plot->plot(xyz, c, press, "xyz-c-press",  time.current_step());
    ns.grad(press);
#endif

    ns.solve(ResRat(1e-3));

#ifdef GHOST
    f = 0.0;
    conc.curvature();
    pr.ghost(c,kappa);
#else
    p=0.0;
#endif

    //if (multigrid.vcycle(ResRat(1e-4))) OMS(converged);
    pr.solve(ResRat(1e-4));
    
#ifdef GHOST
    ns.project_ghost(press,c,kappa);
#else
    ns.project(p);
    press += p;
#endif

    ns.cfl_max();

#if 1
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
    std::cout<<pexact<<"\n";
    for_vijk(press,i,j,k)
      if(c[i][j][k]>=0.99){
        psum += press[i][j][k];
        icount ++;
        perr=pow((press[i][j][k]-pexact),2.0);
      }

    perr=perr/(real(icount)*pow(pexact,2.0));
    perr=sqrt(perr);
    std::cout<<"pressure= "<<psum/real(icount)<<","<<perr<<","<<icount<<","<<"\n";
#endif

    if(time.current_step() % (nint)==0 ||time.current_step()==1 ) {
      boil::plot->plot(uvw, c, press, "uvw-c-press",  time.current_step());
      boil::plot->plot(xyz, c, press, "xyz-c-press",  time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.2 2017/08/15 10:02:45 sato Exp $'/
+-----------------------------------------------------------------------------*/
