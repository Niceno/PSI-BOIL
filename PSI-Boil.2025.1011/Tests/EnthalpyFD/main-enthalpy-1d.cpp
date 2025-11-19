#include "Include/psi-boil.h"
#include <fstream>

/* domain dimensions (given by problem) */
const real LX =   0.5;
const real LY =   0.2;
const real LZ =   0.2;

/* computed parameters */
const int NX = 20;
const int NY = 4;
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
  Grid1D gx( Range<real>( 0,LX), NX, Periodic::yes() );
  Grid1D gy( Range<real>( 0,LY), NY, Periodic::yes());
  Grid1D gz( Range<real>( 0,LZ), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 10;
  const int  nint = 1;
  const real dt  = 0.25 * LX / real(NX);
  Times time(ndt, dt); 
	
  OPR(  NX );
  OPR(  dt );
  OPR( ndt );

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // vel
  Scalar tpr(d), q  (d); // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }
  
  tpr.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  tpr = tpr.shape();
  q = tpr.shape();

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

  for_vijk(tpr,i,j,k){
    if(tpr.xc(i)>0.1 && tpr.xc(i)<0.3){
      tpr[i][j][k] = 1.0;
    } else {
      tpr[i][j][k] = 0.0;
    }
  }

  tpr.exchange_all();

  boil::plot->plot(uvw,tpr, "uvw-tpr", 0);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter water(d);
  water.mu    (1.0000e-3);
  water.rho   (1000.0);
  water.cp    (1.0);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  Enthalpy enth(tpr, q, uvw, time, solver, &water);
  enth.convection_set(TimeScheme::forward_euler());

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "#                       " << boil::endl;
    boil::oout << "########################" << boil::endl;

    enth.new_time_step();
    enth.discretize();
    enth.solve(ResRat(1e-4));

    if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw,tpr, "uvw-tpr",  time.current_step());
    }

  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
