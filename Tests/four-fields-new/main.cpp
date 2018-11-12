/*-------------------------------------------------------------+
| seven bubbles and two droplets merging with a free surface   |
| FSL is used for bubbles and droplets, CIP for free surface   |
+-------------------------------------------------------------*/
#include "Include/psi-boil.h"
#include <iomanip>

void update_step(const Scalar & c, Scalar & step, Scalar & sflag);

const int NX= 48;
const int NZ= NX*2;
const int nint = 200;
const real DB_avg =  0.0015;
const real LX = DB_avg * 6.0;
const real LZ = LX * NZ/NX; 
const real gravity = 9.8;

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
 
  Grid1D gz( Range<real>(0.0, LZ), 
             Range<real>( LZ/NZ,  LZ/NZ ),
              NZ, Periodic::no());

  Grid1D gx( Range<real>(-LX/2.0, LX/2.0),
             Range<real>( LX/NX,  LX/NX ),
                      NX, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter air(d), water(d);
  air.rho   (1.205);
  air.mu    (1.82e-5);
  water.rho   (998.2);
  water.mu    (0.001); 

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());  // ic2//

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); /* velocity */
  Scalar p  (d), f  (d), press(d); /* pressure */

  Scalar c_sep (d, "separated"), g  (d); /* separated phase */
  Scalar c_bub (d, "bubbles"); /* bubbles  */
  Scalar c_dro(d, "droplets"); /* dropltes */

  Scalar step(d), sflag(d); /* step version of c_sep */

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
  
  c_sep = p.shape();
  step  = p.shape();
  sflag = p.shape();
  c_bub = p.shape();
  c_dro = p.shape();

  Matter mixed(water, air, &step, &c_dro, &c_bub);
  mixed.sigma(0.073);

  /*------------+
  |  time step  |
  +------------*/
  const real dxmin = std::min(LX/NX,LZ/NZ);
  const real dt  = 5.0 * pow(0.5*pow(dxmin,3.0)
                  /(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int ndt = 10000; 
  Times time(ndt, dt);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns( uvw, xyz, time, solver, &mixed);
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());

  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(false);
  multigrid.max_cycles(15);

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  for_vijk(c_sep,i,j,k)
    c_sep[i][j][k] = 1.0;
  for_vk(c_sep,k)
    if(c_sep.zc(k) < 0.625 * LZ){ 
      for_vij(c_sep,i,j)
        c_sep[i][j][k]=1.0;
    } else {
      for_vij(c_sep,i,j)
        c_sep[i][j][k]=0.0;
    }

  c_sep.exchange_all();

  update_step(c_sep, step, sflag);  /* convert c_sep to step-function */

  Dispersed bubbles (c_bub, & c_sep, 0, uvw, time, &mixed); 
  Dispersed droplets(c_dro, & c_sep, 1, uvw, time, &mixed); 

  const real r = 0.002; 
  const real a = boil::pi/4.0 + boil::pi/2.0;
  const real DB2 = 0.002;
  const real DB1 = 0.001;

  /* add seven bubbles of diameter 1 and 2 mm */
  bubbles.add( Particle( Position( r * cos(a), 
                             0.125*r * sin(a), 
                                   0.4 * LZ), 
                         Diameter(DB2)));  

  bubbles.add( Particle( Position( -r  * cos(a), 
                                   r  * sin(a), 
                                   0.325 * LZ), 
                         Diameter(DB1), 
                         Position(0.0,0.0,0.1) ) );

  bubbles.add( Particle( Position( -r * cos(a), 
                                   -r * sin(a),
                                    0.17 * LZ), 
                              Diameter(DB2) ) );

  bubbles.add( Particle( Position( r * cos(a), 
                                  -r * sin(a), 
                                   0.5 * LZ), 
                            Diameter(DB1) ) );

  bubbles.add( Particle( Position( -r * cos(a), 
                                0.3*r * sin(a), 
                                   0.45 * LZ), 
                            Diameter(DB1) ) );

  bubbles.add( Particle( Position( -r * cos(a), 
                                   -r * sin(a), 
                                     0.3 * LZ), 
                         Diameter(DB2), 
                         Position(0.08,0.0,0.07) ) );

  bubbles.add( Particle( Position( -1.05 *r * cos(a), 
                                    r * sin(a), 
                                    0.05 * LZ), 
                            Diameter(DB1) ) );

  /* add two droplets of diameter 1 and 2 mm */
  droplets.add( Particle( Position( r * cos(a), 
                                    r * sin(a), 
                                   0.775 * LZ), 
                           Diameter(DB2)));  

  droplets.add( Particle( Position(-r * cos(a), 
                                    r * sin(a), 
                                   0.87 * LZ), 
                           Diameter(DB1)));  

  CIPCSL2 conc (c_sep, g, uvw, time, solver);
  conc.set_nredist(1);
  conc.set_itsharpen(4);
  conc.front_minmax();

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

    bubbles .advance(& xyz);
    droplets .advance(& xyz);

    /*------------------------------------+
    | this is required for CIP since it   |
    | may be modified in Dispersed class  |
    +------------------------------------*/
    c_sep.bnd_update();
    c_sep.exchange_all();
    conc.update_node(c_sep);
 
    /* advance CIP */
    conc.advance();
    conc.totalvol();
    conc.front_minmax();
    update_step(c_sep, step, sflag);
    conc.tension(&xyz, mixed, step);

    /* gravity force */
    const Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] -= gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);

    ns.cfl_max();

    /* essential for moving boundary */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-4));

    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-3));
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();

    /* dt control */
    time.control_dt(ns.cfl_max(),0.25,dt);

    /* output data */
    if(time.current_step() % nint == 0) 
      boil::plot->plot(uvw, c_sep, c_bub, c_dro, "uvw-c", time.current_step());

    /*---------+
    |  exit ?  |
    +---------*/
    if(conc.get_zmaxft()>=LZ*0.99){
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
/***************************************************************************//**
* convert color function to step function
* color function             : c
* step function              : step
* temporary used for flagging: sflag
*******************************************************************************/

  const real phisurf=0.5;
  for_avijk(sflag,i,j,k) {
    sflag[i][j][k]=0;
  }
  for_vijk(c,i,j,k) {
    if(c[i][j][k]>=phisurf)
      sflag[i][j][k]=1;
  }

  // i-direction //
  for(int i=c.si()-1; i<=c.ei(); i++){
    for_vjk(c,j,k){
       if((c[i][j][k]-phisurf)*(c[i+1][j][k]-phisurf)<=0.0){
          sflag[i  ][j][k]=2;
          sflag[i+1][j][k]=2;
       }
    }
  }

  // j-direction //
  for(int j=c.sj()-1; j<=c.ej(); j++){
    for_vik(c,i,k){
      if((c[i][j][k]-phisurf)*(c[i][j+1][k]-phisurf)<=0.0){
          sflag[i][j  ][k]=2;
          sflag[i][j+1][k]=2;
       }
    }
  }

  // k-direction /
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

  step.exchange_all();

}
