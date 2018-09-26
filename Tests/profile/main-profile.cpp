#include "Include/psi-boil.h"

#include <vector>

/* boundary conditions */
const real LX  =   1.0;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), 50, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gx);

  /*-------------------------------+
  |  create a sample scalar field  |
  +-------------------------------*/
  Scalar f(d);

  /*----------------------+
  |  create two profiles  |
  +----------------------*/
  Profile p("profile.dat");
  Profile q(p);
 
  q.scale_values(0.5);
  boil::oout << q << boil::endl;

  /*--------------------------------+
  |  assign values to scalar field  |
  +--------------------------------*/
  for_vijk(f,i,j,k) {
    real r = sqrt( f.xc(i)*f.xc(i) + f.yc(j)*f.yc(j) + f.zc(k)*f.zc(k) );
    f[i][j][k] = q.value_at(r);
  }

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(f, "test", 0);
}	

/*-----------------------------------------------------------------------------+
 '$Id: main-profile.cpp,v 1.1 2009/05/07 17:56:44 niceno Exp $'/
+-----------------------------------------------------------------------------*/
