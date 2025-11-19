#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  boil::aout << "Inside the PSI-Boil!" << boil::endl;

  boil::timer.stop();
}
