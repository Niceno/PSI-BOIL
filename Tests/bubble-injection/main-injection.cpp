#include "Include/psi-boil.h"

#include <vector>
//#include <cmath>
//#include <cstdlib>
//#include <limits>
//#include <ctime>

/* parameters */
const int NX = 60;
const int NY = NX;
const int NZ = 4*NX;
//const int NZ = NX;
const int NINJ = 12;

const real LX = 0.02;
const real LY = LX/real(NX)*real(NY);
const real LZ = LX/real(NX)*real(NZ);

const real gravity=9.8;

/****************************************************************************/
main(int argc, char * argv[]) {

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  boil::timer.start();

  /*--------+
  |  grids  |
  +--------*/
  Grid1D gx(Range<real>( -LX/2.0,      LX/2.0 ), NX, Periodic::no());
  Grid1D gy(Range<real>( -LY/2.0,      LY/2.0 ), NY, Periodic::no());
  Grid1D gz(Range<real>( 0.0,  LZ ), NZ, Periodic::no());

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
  Vector uvw(d), xyz(d); // velocity
  Scalar press(d), p  (d), f  (d); // pressure
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd(Range<int>(NX/2+1-NINJ/2,NX/2+NINJ/2)
                        , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                        , Dir::kmin()
                        , BndType::inlet(),0.0,0.0,0.1 ) );
    uvw.bc(m).add( BndCnd(Range<int>(1,NX)
                        , Range<int>(1,NY/2+1-NINJ/2-1)
                        , Dir::kmin()
                        , BndType::wall() ) );
    uvw.bc(m).add( BndCnd(Range<int>(1,NX/2+1-NINJ/2-1)
                        , Range<int>(NY/2-NINJ/2,NY/2+NINJ/2)
                        , Dir::kmin()
                        , BndType::wall() ) );
    uvw.bc(m).add( BndCnd(Range<int>(NX/2+NINJ/2+1,NX)
                        , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                        , Dir::kmin()
                        , BndType::wall() ) ); 
    uvw.bc(m).add( BndCnd(Range<int>(1,NX)
                        , Range<int>(NY/2+NINJ/2+1,NY)
                        , Dir::kmin()
                        , BndType::wall() ) ); 

    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
  }

  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press = p.shape();

  std::cout<<NX/2-NINJ/2<<" "<<NX/2+NINJ/2<<"\n";
  c.bc().add( BndCnd(Range<int>(NX/2+1-NINJ/2,NX/2+NINJ/2)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd(Range<int>(1,NX)
                   , Range<int>(1,NY/2+1-NINJ/2-1)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd(Range<int>(1,NX/2+1-NINJ/2-1)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd(Range<int>(NX/2+NINJ/2+1,NX)
                   , Range<int>(NY/2+1-NINJ/2,NY/2+NINJ/2)
                   , Dir::kmin()
                   , BndType::wall() ) );
  c.bc().add( BndCnd(Range<int>(1,NX)
                   , Range<int>(NY/2+NINJ/2+1,NY)
                   , Dir::kmin()
                   , BndType::wall() ) );

  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (   1.0e-5  );
  air  .rho   (   1.25    );
  water.mu    (   1.0e-3  );
  water.rho   ( 1000.0    );
  Matter mixed(water,air, &c);
  mixed.sigma(0.072);

  /* Time */
  const real dxmin = std::min(LX/NX,LZ/NZ);
  const real dt  = 5.0*pow(0.5*pow(dxmin,3.0)/(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int nint=100;

  Times time(50000, dt); /* ndt, dt */
  /*---------+
  |  circle  |
  +---------*/
  for_avijk(c,i,j,k)
    c[i][j][k]=1.0;

  Pressure pr( p,   f,   uvw, time, solver, &mixed );
  Momentum ns( uvw, xyz,      time, solver, &mixed );
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::upwind()); 

  CIPCSL2 conc(c, g, kappa, uvw, time, solver);
  conc.set_itsharpen(10);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_cangle(0.0);

  boil::plot->plot(uvw,c,press,"uvw-c-press", 0);

  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    for_m(m)
      for_vmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /* advance */
    conc.new_time_step();
    //conc.convection();
    conc.advance();
    //conc.sharpen();
    conc.totalvol();
    conc.tension(&xyz, mixed);

    ns.discretize();
    pr.discretize();

    /* gravity */
    Comp m=Comp::w();
#if 0
    for(int i=1; i<xyz.ni(m)-1; i++)
      for(int j=1; j<xyz.nj(m)-1; j++)
        for(int k=1; k<xyz.nk(m)-1; k++)
          xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
#else
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] -= gravity * xyz.dV(m,i,j,k)
                       * (mixed.rho(m,i,j,k)-water.rho()->value());
#endif

    /* solve momentum */
    ns.cfl_max();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-6));

    /* solve pressure */
    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-3))) OMS(converged);
    p.exchange();
    ns.project(p);
    press += p;
     press.exchange();
  
    if( time.current_step() % nint == 0 ) {
      boil::plot->plot(uvw,c,press,"uvw-c-press", time.current_step());
    }
//    if( time.current_step() % 500 == 0 ) {
//      uvw.save("uvw", time.current_step());
//      c.  save("c",   time.current_step());
//    }

#if 0
  std::ofstream fout;
  fout.open("profile.txt");
  for_vi(c,i) {
       fout << c.xc(i) <<" "<<uvw[Comp::w()][i][1][1]<<"\n";
  }
  fout.close();
#endif

  }

  boil::timer.stop();
  boil::timer.report();
}
