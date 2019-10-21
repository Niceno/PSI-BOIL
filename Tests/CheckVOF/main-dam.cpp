#include "Include/psi-boil.h"
#define USE_VOF

/* parameters you should set */
const int  LEVEL          = 0;     /* 0,1,2,3 (64,128,256,512 cells) */

/* domain dimensions (given by problem) */
const real inch = 0.254;
const real LX =   9.0*inch;
const real LY =   0.5*inch;
const real LZ =   9.0*inch;

/* computed parameters */
const int NX = 64 * (int)pow(2.0,LEVEL);
const int NZ = NX;

const real gravity = 9.8;

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
  Grid1D gx( Range<real>( 0,LX), NX, Periodic::no() );
  Grid1D gy( Range<real>( 0,LY),  4, Periodic::yes());
  Grid1D gz( Range<real>( 0,LZ), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air  .mu    (1.862e-5);
  air  .rho   (1.1763e+0);
  water.mu    (8.544e-4);
  water.rho   (9.9662e+2);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int ndt = 4000;
  const real dt  = 0.5 * 0.001 / pow(2.0,LEVEL);
  const real tint = 0.1;
  const real tend = 1.0;
  Times time(ndt, dt); 
	
  OPR(  NX );
  OPR(  dt );
  OPR( ndt );

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d), kappa(d); // concentration
  Scalar press(d);

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press = p.shape();
  
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );
  g = c.shape();
  kappa = c.shape();

  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

  const real radius=0.0*inch;
  const real xcent =2.25*inch-radius;
  const real zcent =4.5 *inch-radius;
  for_vijk(c,i,j,k) {
    real dist=pow(c.xc(i)-xcent,2.0)+pow(c.zc(k)-zcent,2.0);
    if (dist<pow(radius*0.75,2.0)) {
      c[i][j][k]=1.0;
      press[i][j][k] = water.rho()->value() * gravity * c.zc(k);
    } else if(dist<pow(radius*1.25,2.0)) {
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
            real dist=pow(xxc-xcent,2.0)+pow(zzc-zcent,2.0);
            if (dist<pow(radius,2.0)){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
      press[i][j][k] = water.rho()->value()*gravity*c.zc(k)*c[i][j][k]
                      +air.rho()->value()*gravity*c.zc(k)*(1.0-c[i][j][k]);
    } else {
      press[i][j][k] = air.rho()->value() * gravity * c.zc(k);
    }
  }

  for_vijk(c,i,j,k) {
    if( (c.xc(i) < xcent && c.zc(k) < 4.5*inch)
      ||(c.xc(i) < 2.25*inch && c.zc(k) < zcent) ){
      c[i][j][k] = 1.0;
      press[i][j][k] = water.rho()->value() * gravity * c.zc(k);
    } else {
      press[i][j][k] = air.rho()->value() * gravity * c.zc(k);
    }
  }
 
  c.exchange();
  boil::plot->plot(uvw, c, press, "uvw-c-press", 0);

  Matter mixed(water, air, &c);
  mixed.sigma (7.6000e-2);

#ifdef USE_VOF
  VOF conc  (c,  g, kappa, uvw, time, solver);

#else
  CIPCSL2 conc  (c,  g, kappa, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
#endif

  conc.front_minmax();
  conc.totalvol();

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /* reset body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* gravity */
    Comp m=Comp::w();
    for(int i=1; i<xyz.ni(m)-1; i++)
      for(int j=1; j<xyz.nj(m)-1; j++)
        for(int k=1; k<xyz.nk(m)-1; k++)
          xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    /* surface tension */
    //conc.tension(&xyz, mixed);
    //boil::plot->plot(xyz, c, press, "xyz-c-press",  time.current_step());

    /* essential for moving front */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(0.001));
    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-6))) OMS(converged);
    ns.project(p);
    press += p;
    press.exchange();

    /* advect color function */
    conc.advance();
    conc.front_minmax();
    conc.totalvol();

    real cflmax = ns.cfl_max();
    boil::oout<<"cfl_max= "<<time.current_time()<<" "<<cflmax<<"\n";
    real cfl_limit=0.25;  // max limit of CFL
    time.control_dt(cflmax, cfl_limit, dt);

    /* tecplot files */
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      boil::plot->plot(uvw,c,press,"uvw-c-press",iint);
    }

    if((time.current_time()) > tend ) {
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }
  }

  boil::timer.stop();
  boil::timer.report();

}
