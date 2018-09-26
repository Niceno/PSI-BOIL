#include "Include/psi-boil.h"

#include <vector>

/******************************************************************************/
main(int argc, char * argv[]) {

  const real b_des = 1.0;
  real       b_new = 0;

  #include "domain.cpp"

  Times time(2000, 0.05); /* ndt, dt */
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // velocity
  Scalar p(d),   f(d);   // pressure
  Scalar q(d);

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  }

  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  
  const Comp m = Comp::u();
  for(int i=xyz.si(m); i<=xyz.ei(m); i++)
  for(int j=xyz.sj(m); j<=xyz.ej(m); j++)
  for(int k=xyz.sk(m); k<=xyz.ek(m); k++)
    xyz[m][i][j][k] = 2.0 * xyz.dV(m,i,j,k);

  Momentum ns(uvw, xyz,      time, solver, & fluid);
  Pressure pr(p,   f,   uvw, time, solver, & fluid);

  AC multigrid( &pr );

  /*------------+
  |  time loop  |
  +------------*/
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
    
    if( multigrid.vcycle(ResRat(0.01)) ) {
      p.exchange();   // is this needed? check!
      ns.project(p);
      uvw.exchange(); // this is not needed. check
      pr.update_rhs();
    }

    OPR(p.min());
    OPR(p.max());

    /*-------+
    |  save  |
    +-------*/
    if( time.current_step() % 500 == 0 ) {
      ns.get_q(&q);
      boil::plot->plot(q,  "q",   time.current_step());
      boil::plot->plot(p,  "p",   time.current_step());
      boil::plot->plot(uvw,"uvw", time.current_step());
      uvw.save("uvw", time.current_step());
      p.save  ("p",   time.current_step());
    }
  }
  ns.get_q(&q);
  boil::plot->plot(q,  "q",   time.current_step()-1);
  boil::plot->plot(p,  "p",   time.current_step()-1);
  boil::plot->plot(uvw,"uvw", time.current_step()-1);
  uvw.save("uvw_final", time.current_step()-1);
  p.save  ("p_final",   time.current_step()-1);
}	

/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.9 2011/05/25 11:09:21 niceno Exp $'/
+-----------------------------------------------------------------------------*/
