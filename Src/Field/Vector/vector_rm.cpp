#include "vector.h"

/******************************************************************************/
void Vector::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* remove file */
  remove(name.c_str());
  
}
