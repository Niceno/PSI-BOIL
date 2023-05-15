#include "Include/psi-boil.h"

//  Coarse grid version of
//  Sato Y, Niceno B. J Comput Phys. 2012;231:3887-95.
//  http://dx.doi.org/10.1016/j.jcp.2012.01.034

/* boundary conditions */
const real LX =   0.03;
const real LZ =   0.03;
const int  NX =   32;
const real radius =   0.01;
const real cangle = 120.0;

/******************************************************************************/
int main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
//boil::plot = new PlotGMV( AsNodes::no(), Buffers::yes() );
  boil::plot = new PlotTEC(  );
//boil::plot = new PlotTEC( AsNodes::yes(), Buffers::yes() );

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 
             NX, Periodic::no());

  Grid1D gz( Range<real>(0,LZ), 
             int(real(NX)*LZ/LX), Periodic::no());

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
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d), kappa(d); // color function

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
  kappa = p.shape();
  
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );
  
  /*---------+
  |  circle  |
  +---------*/
  for_vijk(c,i,j,k)
    c[i][j][k] = 0.0;

  const real zcent =0.0;
  //const real zcent =radius*0.5;
  for_vijk(c,i,j,k) {
    real pi=acos(-1.0);
    real xx=c.xc(i);
    real yy=c.yc(j);
    real zz=c.zc(k);
    real dd = -sqrt(xx*xx+yy*yy+(zz-zcent)*(zz-zcent))+radius;
    real eps = LX/NX*1.5;
    if(dd<-eps){
      c[i][j][k]=0.0;
    } else if (dd<eps){
      c[i][j][k]=0.5+dd/(2.0*eps)+1/(2*pi)*sin(pi*dd/eps);
    } else {
      c[i][j][k]=1.0;
    }
  }
  c.exchange_all();
  boil::plot->plot(uvw, c, press,"uvw-c-press", 0);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (1.0e-3);
  air  .rho   (1.);
  water.mu    (1.0e-3);
  water.rho   (1.0);

  Matter mixed(water, air, & c);
  mixed.sigma (0.01);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const real dxmin = LX/NX;
  //const real dt  = 5.0*pow(air.rho()->value()*pow(dxmin,3.0)  // CIPCSL2 accept 5.0
  //               / (2.0*boil::pi*mixed.sigma()->value()),0.5);// but VOF diverges
  const real dt  = 1.0*pow(air.rho()->value()*pow(dxmin,3.0)    // for VOF 1.0 is used
                 / (2.0*boil::pi*mixed.sigma()->value()),0.5);

  //const int ndt  = 10000;
  //const int nint = 1000;
  const int ndt  = 500;
  const int nint = 50;
  Times time(ndt, dt); /* ndt, dt */
  time.print_time(false);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.min_cycles(3);

#if 0
  CIPCSL2 conc (c,  g, kappa, uvw, time, solver);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_itsharpen(10);
#endif
#if 1
  VOF conc (c,  g, kappa, uvw, time, solver);
  conc.set_wall_curv_method(CurvMethod::DivNorm());
#endif
  conc.set_cangle(cangle);

  /*----------+
  |  gravity  |
  +----------*/
  const Comp m = Comp::w();

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    /* advance */
    conc.advance();
    conc.front_minmax();
    conc.totalvol();

    /* body force */
    for_m(mc){
      for_vmijk(xyz,mc,i,j,k){
        xyz[mc][i][j][k] = 0.0;
    }}
    conc.tension(&xyz, mixed);

    //for_vmijk(xyz,m,i,j,k){
    //  xyz[m][i][j][k] -= 9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
    //}
 
    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.cfl_max();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(0.0001));

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-4));
    p.exchange();
    ns.project(p);
    for_vijk(p,i,j,k){
      press[i][j][k] += p[i][j][k];
    }
    press.exchange();

#if 1
    /* shift pressure */
    real dltp = 0.0;
    if(boil::cart.iam()==0) dltp = press[1][1][1];
    boil::cart.sum_real(&dltp);
    for_vijk(press,i,j,k)
      press[i][j][k] -= dltp;

    int icount=0;
    real psum=0.0;
    for_vijk(press,i,j,k){
      if(c[i][j][k]>=0.99){
        psum += press[i][j][k];
        icount ++;
      }
    }
    boil::cart.sum_int( & icount);
    boil::cart.sum_real( & psum);

    boil::oout<<"x-min= "<< time.current_time()<<" "<<conc.get_xminft()<<" "
              <<conc.get_xmaxft()<<" "<<conc.get_zmaxft()<<" "
              <<psum/real(icount)<<" "<<icount<<"\n";
#endif

    /*------------+
    |  save/plot  |
    +------------*/
    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, press,"uvw-c-press", time.current_step());
      //conc.save("conc", time.current_step());
      //uvw.save ("uvw",  time.current_step());
      //press.save("press", time.current_step());
    }

    /*---------+
    |  exit ?  |
    +---------*/
    std::ifstream infile;
    infile.open ("stop.now", std::ifstream::in);
    if(infile.good()) exit(0);
    infile.close();
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-drop.cpp,v 1.16 2010/03/25 08:15:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
