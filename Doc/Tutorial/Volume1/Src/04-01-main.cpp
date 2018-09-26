#include "Include/psi-boil.h"

const int N = 256;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  real a;

  boil::oout << "Running ..." << boil::endl;

  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++) {
        a = acos(3.14);
      }

  boil::timer.stop();
  boil::timer.report();
}	

/*-----------------------------------------------------------------------------+
 '$Id: 04-01-main.cpp,v 1.3 2008/11/17 19:23:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
