#include "scalar.h"

/******************************************************************************/
void Scalar::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* remove a file */
  remove(name.c_str());
  
}
