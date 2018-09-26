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
  Scalar p(d); // pressure

  /* assign initial values to pressure */
  for_vijk(p,i,j,k)
    p[i][j][k] = (1.0/16.0 )             * 
                 (2.0 + cos(2.0*p.zc(k))) * 
                 (cos(2.0*p.xc(i)) + cos(2.0*p.yc(j)));

  /* plot scalar */
  boil::plot->plot(p, "pressure", 0);

  boil::timer.stop();
  boil::timer.report();
}
	
/*-----------------------------------------------------------------------------+
 '$Id: 06-02-main.cpp,v 1.3 2008/11/17 19:23:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
