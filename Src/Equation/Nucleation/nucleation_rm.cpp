#include "nucleation.h"

/******************************************************************************/
void Nucleation::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  remove(name.c_str());
  
}
