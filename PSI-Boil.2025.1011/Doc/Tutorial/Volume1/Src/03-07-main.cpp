#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  int a;

  APR(a);

  boil::timer.stop();
}
