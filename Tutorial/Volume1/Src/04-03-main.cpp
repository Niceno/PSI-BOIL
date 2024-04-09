#include "Include/psi-boil.h"

const int N = 256;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  real a;
  real b;

  boil::oout << "Running acos ..." << boil::endl;

  /* algorithm using "acos" */
  boil::timer.start("acos algorithm");
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++) {
        a = acos(3.14);
      }
  boil::timer.stop("acos algorithm");

  boil::oout << "Running cos ..." << boil::endl;

  /* algorithm using "cos" */
  boil::timer.start("cos algorithm");
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      for(int k=0; k<N; k++) {
        a = cos(3.14);
      }
  boil::timer.stop("cos algorithm");

  boil::timer.stop();
  boil::timer.report();
}
