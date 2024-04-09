#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  boil::aout << boil::pid << "Hello World!" << boil::endl;

  boil::timer.stop();
}
