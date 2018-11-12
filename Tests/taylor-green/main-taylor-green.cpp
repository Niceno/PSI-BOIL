#include "Include/psi-boil.h"

/* boundary conditions */
const real LX =   6.2831853071796;
const real LY =   LX;
const real LZ =   LX;
  const real mu = 0.01;       // Re =  100
//const real mu = 0.005;      // Re =  200
//const real mu = 0.0025;     // Re =  400
//const real mu = 0.00125;    // Re =  800
//const real mu = 0.000625;   // Re = 1600
//const real mu = 1.0/3000.0; // Re = 3000

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  // boil::plot = NULL; // no plotting
  boil::plot = new PlotGMV();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D g( Range<real>(0,LX),
            Range<real>( LX/128.0, LX/128.0 ),
            128,
            Periodic::yes() ); 

  /*---------+
  |  domain  |
  +---------*/
  Domain d(g, g, g);
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
  fluid.mu(mu);

  Times time(30, 0.01); /* ndt, dt */ /* was 30 */
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p  (d), f  (d); // pressure
  Scalar eps(d);         // dissipation

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr(p,   f,   uvw, time, solver, &fluid);
  Momentum ns(uvw, xyz,      time, solver, &fluid);

  AC multigrid( &pr );

//  uvw.load("uvw");
  for(int i=0; i<p.ni(); i++) {
    for(int j=0; j<p.nj(); j++)
      for(int k=0; k<p.nk(); k++) {
        p[i][j][k] = (1.0/16.0 )             * 
                     (2.0 + cos(2.0*p.zc(k))) * 
                     (cos(2.0*p.xc(i)) + cos(2.0*p.yc(j)));
      }
  }

  const Comp u=Comp::u();
  for(int i=0; i<uvw.ni(u); i++) {
    for(int j=0; j<uvw.nj(u); j++)
      for(int k=0; k<uvw.nk(u); k++) 
        uvw[u][i][j][k] =  sin(uvw.xc(u,i)) * cos(uvw.yc(u,j)) * cos(uvw.zc(u,k));
  }

  const Comp v=Comp::v();
  for(int i=0; i<uvw.ni(v); i++)
    for(int j=0; j<uvw.nj(v); j++)
      for(int k=0; k<uvw.nk(v); k++) 
        uvw[v][i][j][k] = -cos(uvw.xc(v,i)) * sin(uvw.yc(v,j)) * cos(uvw.zc(v,k));

  int it;

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    ns.cfl_max();
    ns.new_time_step();

    ns.solve(ResRat(0.0001));

    for(int i=0; i<p.ni(); i++)
      for(int j=0; j<p.nj(); j++)
        for(int k=0; k<p.nk(); k++)
          p[i][j][k] = 0.0;
    
    multigrid.vcycle(ResRat(0.0001));
    p.exchange();
    ns.project(p);
    uvw.exchange();

/*
    if(it % 50 == 0) {
      uvw.save("uvw", it);
      uvw.plot_par("uvw", it);
      uvw.plot_gmv("uvw", it);
      p.plot_par("p", it);
      f.plot_par("f", it);
    }
*/    
    /*------------+
    |  compute k  |
    +------------*/
    {
      real        dv     = 0.0;
      real        vol    = 0.0;
      real        ke     = 0.0; 
      real        ke_int = 0.0; 
      static real ke_old = 0.0; 

      for_m(m) {
        int li, lj, lk;
        li = lj = lk = 1;

        if(m==Comp::u()) li++;
        if(m==Comp::v()) lj++;
        if(m==Comp::w()) lk++;
      
        for(int i=1; i<uvw.ni(m)-li; i++)
          for(int j=1; j<uvw.nj(m)-lj; j++)
            for(int k=1; k<uvw.nk(m)-lk; k++) {
              dv     =  d.dxc(i) * d.dyc(j) * d.dzc(k);
              ke_int += 0.5 * uvw[m][i][j][k] * uvw[m][i][j][k] * dv;
              vol    += dv;
            }
        } // m
      vol /= 3.0;

      boil::cart.sum_real( &vol );
      boil::cart.sum_real( &ke_int );

      ke = ke_int/vol;

      boil::oout << " K    = " << ke     << boil::endl;
      boil::oout << " K o  = " << ke_old << boil::endl;
      boil::oout << " vol  = " << vol << boil::endl;
      boil::oout << " eps1 = " << (ke_old-ke) / time.dt() << boil::endl;
      ke_old = ke;
    }

    /*--------------+
    |  compute eps  |
    +--------------*/
    {
      ns.get_eps(&eps);       

      real dv      = 0.0;
      real vol     = 0.0;
      real eps_int = 0.0; 

      for(int i=1; i<d.ni()-1; i++)
      for(int j=1; j<d.nj()-1; j++)
      for(int k=1; k<d.nk()-1; k++) {

        dv      =   d.dxc(i) * d.dyc(j) * d.dzc(k);
        eps_int += eps[i][j][k] * dv;
        vol     += dv;

      }
      boil::cart.sum_real( &vol );
      boil::cart.sum_real( &eps_int );

      boil::oout << " vol  = " << vol         << boil::endl;
      boil::oout << " eps2 = " << eps_int/vol << boil::endl;
    }

  } // it

//it--;
//uvw.save("uvw", it);
//uvw.plot_par("uvw", it);
//uvw.plot_gmv("uvw", it);
//p.plot_par("p", it);
//f.plot_par("f", it);

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(p, "test", 0);
}
