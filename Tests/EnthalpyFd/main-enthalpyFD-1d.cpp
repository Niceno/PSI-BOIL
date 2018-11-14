#include "Include/psi-boil.h"
#include <fstream>

/* computed parameters */
const int NX = 100;
const int NY = 4;
const int NZ = 4;

/* domain dimensions (given by problem) */
const real LX =   0.5;
const real LY =   LX*real(NY)/real(NX);
const real LZ =   LX*real(NY)/real(NX);


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
  Grid1D gx( Range<real>( 0,LX)
           , Range<real>(LX/NX/2,LX/NX/2)
           , NX, Periodic::yes() );
  Grid1D gy( Range<real>( 0,LY), NY, Periodic::no());
  Grid1D gz( Range<real>( 0,LZ), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt =  500;
  const int  nint = 50;
  const real dt  = 0.25 * LX / real(NX);
  Times time(ndt, dt); 
	
  OPR(  NX );
  OPR(  dt );
  OPR( ndt );

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // vel
  Scalar c  (d), g  (d); // concentration
  Scalar tpr(d), q  (d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
  }
  
  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );

  tpr.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (1.0000e-5);
  air  .rho   (1.1768e+0);
  air  .cp    (1.0e+3);
  air  .lambda(1.0e-3);
  water.mu    (1.0000e-3);
  water.rho   (1.0000e+3);
  water.cp    (1.0e+3);
  water.lambda(1.0e-3);
  const real latent=1.0e+6;
  const real tsat  =100.00;

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

  uvw.exchange(); // set periodic boundary condition

  for_vijk(c,i,j,k){
    c[i][j][k] = 0.0;
    tpr[i][j][k] = tsat;
  }

  for_vijk(c,i,j,k) {
    if( 0.10 < c.xc(i) && c.xc(i) < 0.20 ){
      c[i][j][k]=1.0;
    tpr[i][j][k] = tsat+1.0;
    }
  }

  
  c.exchange_all();
  tpr.exchange_all();
  boil::plot->plot(uvw,c,tpr, "uvw-c-tpr", 0);

  Matter mixed(water, air, c);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  CIPCSL2 conc (c,  g, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
  //conc.set_epss(0.5);
  conc.front_minmax();
  conc.totalvol();
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, &mixed, tsat);
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

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

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthalpy");

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, tpr, "uvw-c-tpr",  time.current_step());
    }

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
