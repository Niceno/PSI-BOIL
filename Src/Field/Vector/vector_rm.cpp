#include "vector.h"

/******************************************************************************/
void Vector::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* remove file */
  remove(name.c_str());
  
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_rm.cpp,v 1.1 2011/09/22 14:56:43 sato Exp $'/
+-----------------------------------------------------------------------------*/
