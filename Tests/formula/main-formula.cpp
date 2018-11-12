#include "Include/psi-boil.h"

#include <vector>

/* boundary conditions */
const real LX  =   64.0;

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
  Grid1D gx( Range<real>(-LX/2,LX/2), 
             Range<real>(LX/64, LX/64),
             64, Periodic::yes() );

  Formula f;   

  /*-------------+
  |  Example 1:  |
  +-------------*/
  f.evaluate("a=3");          // f.set("a",3); would do the same
  f.evaluate("b=4");          // f.set("b",4); would do the same
  boil::oout << f.evaluate("a+b") << boil::endl;     // prints 7

  /*-------------+
  |  Example 2:  |
  +-------------*/
  Domain d(gx, gx, gx);
  Scalar u(d); 

  /* equation we would like to evaluate */
  std::string equation("x*x + y*y + z*z");

  /* browse through all the cells of variable "u" */
  for_vijk(u, i,j,k) {

    /* set equation variables */
    f.set("x", u.xc(i));
    f.set("y", u.yc(j));
    f.set("z", u.zc(k));

    /* evaluate equation */
    u[i][j][k] = f.evaluate(equation);
  }

  /* used for testing only */
  boil::plot->plot(u, "equation", 0);
}
