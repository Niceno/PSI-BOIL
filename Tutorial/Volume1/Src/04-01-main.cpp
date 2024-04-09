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
        a = acos(-1.0);
      }

  boil::timer.stop();
  boil::timer.report();
}
