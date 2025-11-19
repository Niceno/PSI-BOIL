#include "Include/psi-boil.h"

/* parameters */
const real LX =   1.0;
const real LY =   0.125;

const int NX = 64;
const int NY =  4;

const real Pr = 0.71;
const real Ra = 1.0e+5;

/******************************************************************************/
int main(int argc, char * argv[]) {
  std::cout<<"Hello world\n";
  boil::timer.start();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>( -0.5*LX, 0.5*LX ), NX, Periodic::no());
  Grid1D gy( Range<real>( 0, LY ),           NY, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar t  (d), g  (d); // t.

  /*--------------------+
  |  monitoring points  |
  +--------------------*/
  for_vijk(p,i,j,k){p[i][j][k]=p.dV(i,j,k);}
  Location m0("m0", d, NX, NY, NX);
  m0.print(p);
  //std::cout<<"Rack r4\n";
  Rack    r4("r4", d, Range<int>(1,NX), Range<int>(1,NY), 1);
  //std::cout<<"r4.save_grid\n";
  r4.save_grid(t,"r4-scalar");
  //std::cout<<"main:volume= "<<p.dV(2+NX,2+NY,2+1)<<"\n";
  //std::cout<<"main:dxc= "<<p.dxc(p.si()+NX-1)<<" "<<p.dyc(p.sj()+NY-1)<<" "<<p.dzc(p.sk()+1-1)<<"\n";
  r4.save(p,"r4-volume",0);
  r4.save_grid(uvw,Comp::u(),"r4-vector-U");
  r4.save_grid(uvw,Comp::v(),"r4-vector-V");
  r4.save_grid(uvw,Comp::w(),"r4-vector-W");
  p=0.0;

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
  
  t.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), +0.5 ) );
  t.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), -0.5 ) );
  t.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  t.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  
  /*---------------------------------------+
  |  physical properties, time and solver  |
  +---------------------------------------*/
  Matter fluid(d);
  fluid.mu( Pr );

  //Times time(4000, 0.000075); /* ndt, dt */
  Times time(10, 0.000075); /* ndt, dt */
	
  Krylov * solver = new CG(d, Prec::di());

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr  ( p,   f,   uvw, time, solver, &fluid);
  Momentum ns  ( uvw, xyz,      time, solver, &fluid);
  Enthalpy enth( t,   g,   uvw, time, solver, &fluid);

  AC multigrid( &pr );

  //Location m0("m0", d, NX/2, NY/2, NX/4);
  //Location m0("m0", d, NX, 1, NX);
  //m0.print(t);
  //Rack    r4("r4", d, Range<int>(1,NX), Range<int>(1,NY), 1);
  //r4.save_grid(t,"r4-scalar");
  //r4.save_grid(uvw,Comp::u(),"r4-vector-U");
  //r4.save_grid(uvw,Comp::v(),"r4-vector-V");
  //r4.save_grid(uvw,Comp::w(),"r4-vector-W");
  //Rack    r5("r5", d, Range<int>(1,NX), 1, 1);
  //std::cout<<"call r5.save_grid\n";
  //r5.save_grid(t,"r5-t");
  //r5.save_grid(uvw,Comp::u(),"r5-U");

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;
	  
    enth.new_time_step();
    enth.solve(ResRat(0.001));

    ns.cfl_max();
    ns.new_time_step();

    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k)
      xyz[m][i][j][k] = Pr*Ra * 0.5*(t[i][j][k]+t[i][j][k-1]) * xyz.dV(m,i,j,k);

    ns.solve(ResRat(0.001));

    multigrid.vcycle(ResRat(0.001));
    p.exchange();
    ns.project(p);

    pr.update_rhs();

    //m0.print(uvw,Comp::u());
    m0.print(t);
  }

  boil::plot = new PlotTEC();
  boil::plot->plot(uvw, t, "uvw-t", time.current_step()-1);
  boil::plot->plot(p,      "p",     time.current_step()-1);

  //Rack    r0("u-comp", d, NX/2+1, NY/2, Range<int>(1,NX));
  //Rack    r2("w-comp", d, Range<int>(1,NX), NY/2, NX/2+1);
  //Rack    r3("r3", d, Range<int>(1,NX), 1, 1);
  //Rack    r4("r4", d, Range<int>(1,NX), Range<int>(1,NY), 1);
  //r0.print(uvw,Comp::u());
  //r2.print(uvw,Comp::w());
  //r3.print(t);
  //r3.print(uvw,Comp::u());
  r4.save(t,"r4-tpr",time.current_step());
  r4.save(uvw,Comp::u(),"r4-U",time.current_step());
  r4.save(uvw,Comp::v(),"r4-V",time.current_step());
  r4.save(uvw,Comp::w(),"r4-W",time.current_step());

  //const int ii = t.si();   /* i inside the domain */
  //const int iw = ii-1;     /* i in the wall */
  //const int j  = t.ej()/2; /* j in the middle */
  //boil::oout << " Nusselt number " << boil::endl;
  //for_vk(t,k) {
  //  real nu = (t[iw][j][k] - t[ii][j][k]) / t.dxw(ii);
  //  boil::oout << t.zc(k) << " " << nu << boil::endl;
  //}

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}
