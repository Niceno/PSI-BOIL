#include "Include/psi-boil.h"

const real L =  6.2831853071796;
const int  N = 64;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /* plot in Tecplot format */
  boil::plot = new PlotTEC();

  /* grid in "x", "y" and "z" direction */
  Grid1D g( Range<real>(0,L), N, Periodic::yes() );

  /* computational domain */
  Domain d(g, g, g);
	
  /* define unknowns */
  Vector uvw(d); // velocity

  /* assign initial values to velocity in "i" direction */
  Comp m = Comp::u();
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k] =  sin(uvw.xc(m,i)) * cos(uvw.yc(m,j)) * cos(uvw.zc(m,k));

  /* assign initial values to velocity in "j" direction */
  m = Comp::v();
  for_vmijk(uvw,m,i,j,k)
    uvw[m][i][j][k] = -cos(uvw.xc(m,i)) * sin(uvw.yc(m,j)) * cos(uvw.zc(m,k));

  /* resolution of "u" velocity component */
  OPR(uvw.ni( Comp::u() ));
  OPR(uvw.nj( Comp::u() ));
  OPR(uvw.nk( Comp::u() ));
  /* resolution of "v" velocity component */
  OPR(uvw.ni( Comp::v() ));
  OPR(uvw.nj( Comp::v() ));
  OPR(uvw.nk( Comp::v() ));
  /* resolution of "w" velocity component */
  OPR(uvw.ni( Comp::w() ));
  OPR(uvw.nj( Comp::w() ));
  OPR(uvw.nk( Comp::w() ));

  /* plot vector */
  boil::plot->plot(uvw, "velocity", 0);

  boil::timer.stop();
  boil::timer.report();
}
