/*-------------------------------------------------------------+
|            Simulation of mini - bubble column                |
|           In order to retrieve Deen bubble column:           |
|              LX should be 0.15 instead of 0.05               |
+-------------------------------------------------------------*/

#include "Include/psi-boil.h"
#define resume false  /* this is for resuming the simulation */

/* boundary conditions */
 const int NX = 48; //152;
 const int NZ = NX * 3;
 const real RB = 0.5 * 0.004; 
 const real DB = 2.0 * RB;
 const real LX = 0.05; //0.15; 
 const real LZ = LX*NZ/NX; 

 const int nholes = 2; // 7;  /* this will give 7 x 7 = 49 spargers */
 const real pitch = 0.00625; /* distance between two spargers */
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
  Matter air(d), water(d);
  air.rho   (1.205);
  water.rho (998.2);
  air.mu    (1.82e-5);   
  water.mu  (1.0e-3); 

 /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d), de(d), uli(d); 
  Scalar p  (d), f  (d); 
  Scalar press(d);
  Scalar c_dro(d); /* color function for the droplets */
  Scalar c  (d); /* color function for the bubbles */
  Scalar c_sep(d, "separated"); /* separated phase */
	
  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ));
  }
  for_m(m) {
    uli.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uli.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uli.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uli.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uli.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uli.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ));
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  press = p.shape();  

  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  c_sep = p.shape();
  c_dro = p.shape();

  real dt = 1.0e-4; 
  //int n = int(120.0/dt);
  int n = int(6.0/dt);
  Times time(n, dt); /* ndt, dt */ 

  /* this is for plotting the data */
  //const int  nsint = 120; /* 120: plot every 1 sec; 480 --> 0.25 seconds */
  const int  nsint = 12; /* 120: plot every 1 sec; 480 --> 0.25 seconds */
        int  csint =   0;
  std::vector<real> save_instants;
  save_instants.resize(nsint+1); /* to store last too */
  for(int i=0; i<=nsint; i++) {
    save_instants[i] = (real)i * time.total_time() / (real)nsint;
  }

 /* this is for saving the data */
  const int         nb_bck  = 240;
        int         bck_ind =   0;
  std::vector<real> bck_instants;
  bck_instants.resize(nb_bck + 1); 
  for(int i = 0; i <= nb_bck; i++) {
    bck_instants[i] = (real)i * time.total_time() / (real)nb_bck;
  }

  Matter mixed(air,water, &c_sep, &c, & c_dro);
  mixed.sigma(0.072);

  c_sep = 0.0; 
  c_sep.exchange_all();
  c = 0.0; 
  c.exchange_all();
  c_dro = 1.0; 
  c_dro.exchange_all();

  Dispersed disp (c, & c_sep, 1, uvw, time, &mixed); 

  const real bubble_period = 0.01489; /* period for bubbles release */
  std::vector<real> bubble_instants;
  const int  nbint = time.total_time()/bubble_period;
        int  cbint = 0;
  bubble_instants.resize(nbint+1); /* to store last too */
  for(int i=0; i<=nbint; i++) {
      bubble_instants[i] = (real)i * bubble_period;
  }

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Momentum ns(uvw, xyz, time, solver, &mixed);  
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(ConvScheme::superbee());
  Pressure pr(p, f, uvw, time, solver, &mixed); 
  AC multigrid( &pr );

  multigrid.stop_if_diverging(false);
  multigrid.max_cycles(30);

  /*Locate a point: (0.0 ; 0.0; 0.56 * LZ); */
  real xl = 0.0; real yl = 0.0; real zl = 0.56 * LZ;
  int ixl = d.local_i( d.I(xl) );
  int jyl = d.local_j( d.J(yl) );
  int kzl = d.local_k( d.K(zl) );
 
 #if resume
/* set resume to true when you want to resume the simulation */
  /* the numbers for csint, cbint, bck_ind: you find them in the file output
    Search for the word: SAVING */
  time.first_step(20525); // the step from the ouput file where you have SAVING
  csint= 7;
  cbint= 403;
  bck_ind= 13;
  time.current_time(6.00041); /* select the time of the next
                                  time step after SAVING */
  uvw.load("uvw",12); /* usually this is equal to bck_ind - 1 */
  press.load("press",12);
  disp.load("particles",12);
  time.set_dt(0.00021);/* this can be the difference between current time of 
                          time step after Saving and the one at the SAVING */  
  #endif
  ///////////////////////////////////////////////////////////////////
  //                     Start Time Loop                           //
  ///////////////////////////////////////////////////////////////////
  for(time.start(); time.end(); time.increase()) {

    for_m(m) for_vmijk(xyz,m,i,j,k) xyz[m][i][j][k] = 0.0;

    /*----------+
    |  bubbles  |
    +----------*/
    if(cbint <= nbint) {
      if(time.current_time() > bubble_instants[cbint]) {
        real xb0 = -0.5 * (nholes - 1) * pitch;
        real yb0 = -0.5 * (nholes - 1)* pitch;
        real zb0 = 0.0025; 
        for (int i = 0; i < nholes; i++) {
          for (int j = 0; j < nholes; j++) {
            real xb = xb0 + real(i) * pitch;
            real yb = yb0 + real(j) * pitch;
            disp.add(Particle( Position(xb, yb, zb0), 
                               Diameter(0.004),
                               Position(0.0,0.0,0.4)));
          }
        }
        cbint++;
        boil::oout<<"ADD BUBBLES.. " << boil::endl;
        disp.check_add_particles();
     }
    }

    disp.advance(& xyz);    

    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] -= 9.81 * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
    }

    ns.discretize();
    pr.discretize();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-4));  
    p = 0.0;
    p.exchange();

    multigrid.vcycle(ResRat(1e-3));
    p.exchange();
    ns.project(p);
    press +=p;
    press.exchange();
    pr.update_rhs();

    time.control_dt(ns.cfl_max(), 0.5, 1.0); 

    if(csint <= nsint) {
      if(time.current_time() > save_instants[csint]) {
        boil::plot->plot(uvw, c, "uvw-c-press",  csint);
        csint++;
      }
    }

  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();
}	
