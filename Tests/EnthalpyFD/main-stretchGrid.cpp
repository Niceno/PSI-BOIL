#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
using namespace std;

#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

/* domain dimensions */
//const real LX = 0.005;
const real LX = 0.04;
const real LZ = 1.0;
const int gLevel = 1;  //grid level=2,3,4
const int NX =  4*gLevel;
const int NZ = 20*gLevel;

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
  Grid1D gx ( Range<real>(-LX*0.5,LX*0.5), NX, Periodic::yes() );
  Grid1D gz ( Range<real>(-LZ,1e-8), 
              Range<real>(0.1,0.01),
              NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Domain d(gx, gx, gz, & floor);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d);           // vel
  Scalar p  (d), f  (d), press(d); // pressure
  Scalar c  (d), g  (d);           // concentration
  Scalar tpr(d), q  (d), step(d); // temperature

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  g = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );
  q = tpr.shape();

  /* heater power */
  real qsrc=0.0;
  real qflux= 100.0*1000.0;     // [W/m2]
  if(approx(tpr.zn(1), -LZ)){
    qsrc=qflux/tpr.dzc(1);  // [W/m3]
  }
  boil::oout<<"#qsrc= "<<qsrc<<"\n";

    q=0.0;
    if(approx(tpr.zn(1), -LZ)){
      for_vij(tpr,i,j){
        q[i][j][1]=qsrc*tpr.dV(i,j,1);
      }
    }
    q.exchange();


  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d), copper(d);
  vapor  .mu    (1.255e-5);
  vapor  .rho   (1000.0);
  vapor  .cp    (10000.0);
  vapor  .lambda(100.0);
  liquid.mu    (0.28e-3);
  liquid.rho   (1000.0);
  liquid.cp    (10000.0);
  liquid.lambda(100.0);
  copper.rho   (1000.0);
  copper.cp    (10000.0);
  copper.lambda(100.0);
  Matter mixed(liquid, vapor, c);
  mixed.sigma(5.9e-2);

  real omass_sumt=0.0;

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 200;
  const int  nint = 10;
  const real dxmin = d.dxyz_min();
  boil::oout<<"main:dxmin= "<<dxmin<<" "<<boil::cart.iam()<<"\n";
  real diff_num = 1.0;   // diffusion number
  real diff_coef = copper.lambda()->value()/copper.cp()->value();
  const real dt = diff_num*dxmin*dxmin/diff_coef;
  cout<<"dt= "<<dt<<"\n";
  Times time(ndt, dt);
  time.print_time(false);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::ic2());

  /* enthalpy equation */
  cout<<"Enthal\n";
  real tt = 0.0;
  EnthalpyFD enthFD(tpr, q, c, uvw, time, solver, & mixed , tt, & copper);
  cout<<"Enthal\n";
  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  cout<<"c init\n";
  c = 1.0;
  c.bnd_update();
  c.exchange_all();

  cout<<"tpr init\n";
  tpr = 0.0;
  tpr.bnd_update();
  tpr.exchange_all();
  boil::plot->plot(tpr, "tpr", 0);

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

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    cout<<"discretize\n";
    enthFD.discretize();
    cout<<"solve\n";
    enthFD.solve(ResRat(1e-16),"enthFD");

    /*--------------+
    |  output data  |
    +--------------*/
    if((time.current_step() % nint) == 0 ) {
      boil::plot->plot(tpr,"tpr", time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();
}
