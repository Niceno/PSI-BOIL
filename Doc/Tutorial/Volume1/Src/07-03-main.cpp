#include "Include/psi-boil.h"

const real L =  1.0;
const int  N = 64; 

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  Grid1D g( Range<real>(0,L), N, Periodic::no());
  Domain d(g, g, g);

  Scalar t(d), q(d);                          /* temperature and its source */
  Vector uvw(d);                              /* velocity field */

  t.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 300.0 ) ); /* b.c. */
  t.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(), 300.0 ) ); /* b.c. */
  t.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );          /* b.c. */
  t.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );          /* b.c. */
  t.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );          /* b.c. */
  t.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );          /* b.c. */

  Matter solid(d);                                  /* matter */
  solid.lambda(0.5);

  Krylov * solver = new CG(d, Prec::di());          /* linear solver */

  Times time;                                       /* simulation time */
	
  Enthalpy enth(t, q, uvw, time, solver, &solid);   /* enthalpy conservation
                                                       equation              */
  enth.diffusion_set(TimeScheme::backward_euler()); /* time stepping scheme */

  AC multigrid( &enth );                            /* AMG solver for enth. */ 

  t = 310.0;                                        /* initial "guess" */

  for_vijk(q,i,j,k)                                 /* fill the source term */
    q[i][j][k] = 100.0*q.dV(i,j,k);

  multigrid.vcycle(ResRat(1e-4));                   /* solve linear system */

  boil::plot = new PlotTEC();
  boil::plot->plot(t, "t");                         

  boil::timer.stop();
  boil::timer.report();
}
